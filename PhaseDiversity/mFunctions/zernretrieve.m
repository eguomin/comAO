function [cEstimate, imgEstimate, J] = zernretrieve(imgs, r, theta, idx, p, c_delta, c0,opidx, gamma, it, GPUflag)
% retrieve the zernike coefficients of aberrated phase based on phase
% diversity images;
% update formula(recurrence): 
%                       c_i+1 = c_i + c_delta
%                       c_delta = -H^(-1)g
% Output
%   cEstimate: retrieved phase coefficients(um)
%   imgEstimate: retrieved image
%   J: cost function values at each iteration
% Input
%   imgs: phase diversity images, third dimension: N, the number of images
%   r,theta,idx: define the pupil aperture of the wavefront
%   1)r: a vector of numbers between 0 and 1
%   2)theta: a vector of angles, has same length with r
%   3)idx:
%       elements should be integers(>=0)
%   p: a vector of single indexes(OSA/ANSI convention) for Zernike components,
%   c_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (should be normalized to phase unit: pi); if matrix,
%   the second dimension should be N-1
%   c0: initial coefficients vector (should be normalized to phase unit: pi)
%   opidx: indexes that indicate which coefficients in p need to be updated
%   gamma: parameter for the regularization
%   it: iteration number
%   GPUflag: GPU options, 0: CPU; 1: GPU;

% By: Min Guo
% July 26, 2017
% Modification(Sep. 11, 2017, Min Guo)
% : to enable more than 2 diversity images and feed back cost function
% values
% Modification(Sep. 13, 2017, Min Guo)
% : feed back reconstruct object image
% Modification(Sep. 19, 2017, Min Guo)
% : update the coefficient selectively
% Modification(Sep. 25, 2017, Min Guo)
% : incorporate the regularization term for phase(quadratic penalty term)
% Modification(Sep. 27, 2017, Min Guo)
% : enable GPU calculation
% Modification(Jan. 28, 2020, Min Guo)
% 1) change function name from rePD_v6 to zernretrieve
% 2) change arguments interface
% 3) make zernPolynomials (the basis components of Zern) as a optional
% argrument

iteration = it;
% gamma = 1e-14;
alpha = 0;

J = zeros(iteration, 1);
[Sx, Sy, imgNum] = size(imgs); % image numbers: 2
% zernNum = length(p);
opNum = length(opidx);
deltaNum = size(c_delta,2);
if (imgNum - deltaNum)~= 1
    error('zernretrieve:NMlength','deltaNum should be: imgNum -1.')
end
idxNum = length(r);

imgDs = zeros(Sx, Sy, imgNum); % FT of imgs
for k = 1:imgNum
    img = squeeze(imgs(:,:,k));
    imgD = fftshift(fft2(ifftshift(img)));
    imgDs(:,:,k) = imgD;
end    

zernPolynomials = create_zernpolybasis(p,r,theta); % Zernike polynomial basis

% Calculate Rc
zernP2D = zeros(Sx,Sy);
zernPhs = zeros(idxNum,opNum);
zernPvs = zeros(idxNum,opNum);
Rc = zeros(opNum, opNum);
for m = 1:opNum
    zernM = opidx(m);
    zernP2D(idx) = zernPolynomials(:,zernM);
    zernP2Dshifted = imtranslate(zernP2D,[1 0]) - zernP2D;
    zernPhs(:,m) = zernP2Dshifted(idx);
    zernP2Dshifted = imtranslate(zernP2D,[0 1]) - zernP2D;
    zernPvs(:,m) = zernP2Dshifted(idx);
end
for m = 1:opNum
    for n = 1:opNum
        Rc(m,n) = dot(zernPhs(:,m),zernPhs(:,n)) + dot(zernPvs(:,m),zernPvs(:,n));
    end
end

pupilMask = zeros(Sx, Sy);
pupilMask(idx) = 1;

waveFront_deltas = zeros(Sx,Sy,imgNum); % known aberration, e.g., defocus
waveFront_delta = zeros(Sx,Sy);
for k = 2:imgNum
    waveFront_delta(idx) = create_wavefront(p, c_delta(:,k-1), r, theta);
    waveFront_deltas(:,:,k) = waveFront_delta;
end

% c0 = zeros(zernNum,1);
Hks = zeros(Sx, Sy, imgNum);
hks = zeros(Sx, Sy, imgNum);
Sks = zeros(Sx, Sy, imgNum);
cEstimate = c0;
opEstimate = c0(opidx);
waveFront = zeros(Sx,Sy);
zernPolynomial2D = zeros(Sx, Sy); % matrix corresponding to image size

g = zeros(opNum,1);
DQks = zeros(Sx, Sy, imgNum);
HGN_phi_ns = zeros(Sx,Sy,opNum);
Hmatrix = zeros(opNum,opNum);
% GPU Data transfering
if(GPUflag==1)
    J = gpuArray(J);
    imgDs = gpuArray(imgDs);
    zernPolynomials = gpuArray(zernPolynomials);
    Rc = gpuArray(Rc);
    pupilMask = gpuArray(pupilMask);
    waveFront_deltas = gpuArray(waveFront_deltas);
    
    Hks = gpuArray(Hks);
    hks = gpuArray(hks);
    Sks = gpuArray(Sks);
    cEstimate = gpuArray(cEstimate);
    opEstimate = gpuArray(opEstimate);
    waveFront = gpuArray(waveFront);
    zernPolynomial2D = gpuArray(zernPolynomial2D);
   
    g = gpuArray(g);
    DQks = gpuArray(DQks);
    HGN_phi_ns = gpuArray(HGN_phi_ns);
    Hmatrix = gpuArray(Hmatrix);
end

for i = 1:iteration
    % ================== %
    % Gradient: g
    % ================== %
    c = cEstimate;
    op = opEstimate;
    waveFront(idx) = create_wavefront(p, c, r, theta);
    for k = 1:imgNum
        phi = waveFront+waveFront_deltas(:,:,k); %phase
        Hk = pupilMask.*exp(1i*phi); % pupil function
        hk = fftshift(ifft2(ifftshift(Hk))); % prh
        sk = abs(hk).^2; % PSF
        Sk = fftshift(fft2(ifftshift(sk))); % FT of PSF
        Hks(:,:,k) = Hk;
        hks(:,:,k) = hk;
        Sks(:,:,k) = Sk;       
    end
    
    % calculate F
    temp1 = conj(Sks).*imgDs;
    numerator = sum(temp1,3);
    temp2 = abs(Sks).^2;
    denominator = sum(temp2,3)+gamma;
    F = numerator./denominator;
    
    % calculate Vk and g[phi]
    g_phi = zeros(idxNum,1);
    if(GPUflag==1)
        g_phi = gpuArray(g_phi);
    end
    for k = 1:imgNum
        Sk = squeeze(Sks(:,:,k));
        Dk = squeeze(imgDs(:,:,k));
        Vk = conj(F).*Dk - abs(F).^2.*Sk;
        Vk_iFT = fftshift(ifft2(ifftshift(Vk)));
        hk = hks(:,:,k);
        Hk = Hks(:,:,k);
        temp1 = hk.*real(Vk_iFT);
        temp2 = fftshift(fft2(ifftshift(temp1)));
        temp3 = imag(conj(Hk).*temp2);
        g_phi = g_phi - 2*temp3(idx);
    end
    
    % calculate gradient: g
    for m = 1:opNum
        zernM = opidx(m);
        zernPolynomial = zernPolynomials(:,zernM);
        g(m) = dot(g_phi,zernPolynomial);
    end
    g = g + 0.5 * alpha * Rc * c;
    
    % ======================== %
    % Hessian matrix: Hmatrix
    % ======================== %
    % simplify to K=2: k=2,j=1
    Q = denominator;
    Qsqrt = Q.^0.5;
    for k = 1:imgNum
        Dk = squeeze(imgDs(:,:,k));
        DQk = Dk./Qsqrt;
        DQks(:,:,k) = DQk;
    end
    for n = 1:opNum
        zernN = opidx(n);
        zernPolynomial2D(idx) = zernPolynomials(:,zernN);
        HGN_phi_n = 0;
        for k = 2:imgNum
            Hk = squeeze(Hks(:,:,k));
            hk = squeeze(hks(:,:,k));
            DQk = squeeze(DQks(:,:,k));
            for j = 1:k-1
                % calculate U
            
                Hj = squeeze(Hks(:,:,j));
                hj = squeeze(hks(:,:,j));
                DQj = squeeze(DQks(:,:,j));
                
                Hk_phi_n = Hk.* zernPolynomial2D;
                temp1 = fftshift(ifft2(ifftshift(Hk_phi_n)));
                temp2 = imag(conj(hk).*temp1);
                temp3 = fftshift(fft2(ifftshift(temp2)));
                
                Hj_phi_n = Hj.* zernPolynomial2D;
                temp4 = fftshift(ifft2(ifftshift(Hj_phi_n)));
                temp5 = imag(conj(hj).*temp4);
                temp6 = fftshift(fft2(ifftshift(temp5)));
                
                Ujk = DQj.*temp3 - DQk.*temp6;
                
                % calculate HGN_phi_n
                temp1 = fftshift(ifft2(ifftshift(conj(DQk).*Ujk)));
                temp2 = fftshift(fft2(ifftshift(hj.*temp1)));
                temp3 = conj(Hj).*temp2;
            
                temp4 = fftshift(ifft2(ifftshift(conj(DQj).*Ujk)));
                temp5 = fftshift(fft2(ifftshift(hk.*temp4)));
                temp6 = conj(Hk).*temp5;
            
                HGN_phi_n = 4* imag(temp3-temp6) + HGN_phi_n;    
            end
        end
        HGN_phi_ns(:,:,n) = HGN_phi_n;
    end
        
    % calculate Hessian Matrix elements
    for m = 1:opNum
        zernM = opidx(m);
        zernPolynomial = zernPolynomials(:,zernM); % vector corresponding to idx
        for n = 1:opNum
            HGN_phi_n = squeeze(HGN_phi_ns(:,:,n));
            HGN_phi_n_idx = HGN_phi_n(idx);
            
            Hmatrix(m,n) = dot(HGN_phi_n_idx, zernPolynomial);
        end
    end
    
    Hmatrix = Hmatrix + alpha*Rc;
    
    opEstimate = op - Hmatrix\g;
    cEstimate(opidx) = opEstimate;
    
    
    temp1 = abs(imgDs).^2;
    temp1 = sum(temp1,3);
    temp2 = conj(imgDs).*Sks;
    temp2 = sum(temp2,3);
    temp2 = abs(temp2).^2./Q;
    temp3 = temp1-temp2;
    J(i) = sum(temp3(:));
    
    imgEstimate = real(fftshift(ifft2(ifftshift(F))));
    imgEstimate(imgEstimate<0) = 0;
end

cEstimate = cEstimate/length2phase; % transform phase to length
if (GPUflag == 1)
    cEstimate = gather(cEstimate);
    imgEstimate = gather(imgEstimate);
    J = gather(J);
end


