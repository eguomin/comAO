function [cEstimate, imgEstimate, J] = zernretrieve_loop(imgDs, r, theta, idx,...
    p, waveFront_deltas, c0, zernPolynomials, Rc, gamma, it, GPUflag)
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
%   imgDs: FT of phase diversity images, third dimension: N, the number of images
%   r,theta,idx: define the pupil aperture of the wavefront
%   1)r: a vector of numbers between 0 and 1
%   2)theta: a vector of angles, has same length with r
%   3)idx:
%       elements should be positive integers(>=1)
%   p: a vector of single indexes(Fringe convention) for Zernike components,
%   waveFront_deltas: wavefronts corresponding to phase diversity images
%   (N)
%   Rc: Rc matrix
%   gamma: parameter for the regularization
%   it: iteration number
%   GPUflag: GPU options, 0: CPU; 1: GPU;
% By: Min Guo
% Jan 29, 2020

alpha = 0;
J = zeros(it, 1);
[Sx, Sy, imgNum] = size(imgDs); % image numbers: 2
opNum = length(p);
idxNum = length(r);


% GPU Data transfering
if(GPUflag==1)
    pupilMask = zeros(Sx, Sy, 'gpuArray');
    waveFront = zeros(Sx,Sy, 'gpuArray');
    Hks = zeros(Sx, Sy, imgNum, 'gpuArray');
    hks = zeros(Sx, Sy, imgNum, 'gpuArray');
    Sks = zeros(Sx, Sy, imgNum, 'gpuArray');
    cEstimate = gpuArray(c0);
    opEstimate = gpuArray(c0);
    zernPolynomial2D = zeros(Sx, Sy, 'gpuArray');
    
    g = zeros(opNum,1, 'gpuArray');
    DQks = zeros(Sx, Sy, imgNum, 'gpuArray');
    HGN_phi_ns = zeros(Sx,Sy,opNum, 'gpuArray');
    Hmatrix = zeros(opNum,opNum, 'gpuArray');
    J = gpuArray(J);
else
    pupilMask = zeros(Sx, Sy);
    waveFront = zeros(Sx,Sy);
    
    Hks = zeros(Sx, Sy, imgNum);
    hks = zeros(Sx, Sy, imgNum);
    Sks = zeros(Sx, Sy, imgNum);
    
    cEstimate = c0;
    opEstimate = c0;
    zernPolynomial2D = zeros(Sx, Sy);
    g = zeros(opNum,1);
    DQks = zeros(Sx, Sy, imgNum);
    HGN_phi_ns = zeros(Sx,Sy,opNum);
    Hmatrix = zeros(opNum,opNum);
    
end
pupilMask(idx) = 1;

for i = 1:it
    % ================== %
    % Gradient: g
    % ================== %
    c = cEstimate;
    op = opEstimate;
    waveFront(idx) = create_wavefront(p, c, r, theta);
    for k = 1:imgNum
        phi = waveFront+waveFront_deltas(:,:,k); %phase
        Hk = pupilMask.*exp(1i*phi); % pupil function
        hk = ifft2(Hk); % prh
        sk = abs(hk).^2; % PSF
        Sk = fft2(sk); % FT of PSF
        Hks(:,:,k) = Hk;
        hks(:,:,k) = hk;
        Sks(:,:,k) = Sk;       
    end
    
    % calculate F
    temp1 = conj(Sks).*imgDs;
    numerator = sum(temp1,3);
    temp2 = abs(Sks).^2;
    denominator = sum(temp2,3) + gamma;
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
        Vk_iFT = ifft2(Vk);
        hk = hks(:,:,k);
        Hk = Hks(:,:,k);
        temp1 = hk.*real(Vk_iFT);
        temp2 = fft2(temp1);
        temp3 = imag(conj(Hk).*temp2);
        g_phi = g_phi - 2*temp3(idx);
    end
    
    % calculate gradient: g
    for m = 1:opNum
        zernPolynomial = zernPolynomials(:,m);
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
        zernPolynomial2D(idx) = zernPolynomials(:,n);
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
                temp1 = ifft2(Hk_phi_n);
                temp2 = imag(conj(hk).*temp1);
                temp3 = fft2(temp2);
                
                Hj_phi_n = Hj.* zernPolynomial2D;
                temp4 = ifft2(Hj_phi_n);
                temp5 = imag(conj(hj).*temp4);
                temp6 = fft2(temp5);
                
                Ujk = DQj.*temp3 - DQk.*temp6;
                
                % calculate HGN_phi_n
                temp1 = ifft2(conj(DQk).*Ujk);
                temp2 = fft2(hj.*temp1);
                temp3 = conj(Hj).*temp2;
            
                temp4 = ifft2(conj(DQj).*Ujk);
                temp5 = fft2(hk.*temp4);
                temp6 = conj(Hk).*temp5;
            
                HGN_phi_n = 4* imag(temp3-temp6) + HGN_phi_n;    
            end
        end
        HGN_phi_ns(:,:,n) = HGN_phi_n;
    end
        
    % calculate Hessian Matrix elements
    for m = 1:opNum
        zernPolynomial = zernPolynomials(:,m); % vector corresponding to idx
        for n = 1:opNum
            HGN_phi_n = squeeze(HGN_phi_ns(:,:,n));
            HGN_phi_n_idx = HGN_phi_n(idx);
            
            Hmatrix(m,n) = dot(HGN_phi_n_idx, zernPolynomial);
        end
    end
    
    Hmatrix = Hmatrix + alpha*Rc;
    
    opEstimate = op - Hmatrix\g;
    cEstimate = opEstimate;
    
    
    temp1 = abs(imgDs).^2;
    temp1 = sum(temp1,3);
    temp2 = conj(imgDs).*Sks;
    temp2 = sum(temp2,3);
    temp2 = abs(temp2).^2./Q;
    temp3 = temp1-temp2;
    J(i) = sum(temp3(:));
    
    imgEstimate = real(ifft2(F));
    imgEstimate(imgEstimate<0) = 0;
end

if (GPUflag == 1)
    J = gather(J);
    cEstimate = gather(cEstimate);
    imgEstimate = gather(imgEstimate);
end


