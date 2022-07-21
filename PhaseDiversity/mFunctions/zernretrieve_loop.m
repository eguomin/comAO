function [cEstimate, imgEstimate, J, itTotal, cInter] = zernretrieve_loop(imgDs, ...
    r0, theta0, idx, p, waveFront_deltas, c0, zernPolynomials, Rc, ...
    penalChioce, gamma, alpha, itLimit, stopChoice, tolValue, GPUflag, nFlag)
% retrieve the zernike coefficients of aberrated phase based on phase
% diversity images;
% update formula(recurrence): 
%                       c_i+1 = c_i + c_delta
%                       c_delta = -H^(-1)g
% Output
%   cEstimate: retrieved Zernike coefficients(um)
%   imgEstimate: retrieved image
%   J: cost function values at each iteration
%   itTotal: actual iteration number
%   cInter: Zernike coefficients at each iteration
% Input
%   imgDs: FT of phase diversity images, third dimension: N, the number of images
%   r0,theta0,idx: define the pupil aperture of the wavefront
%   1)r0: a 2D matrix of numbers between 0 and 1
%   2)theta0: a 2D matrix of angles (rad), has same size with r0
%   3)idx:
%       elements should be positive integers(>=1)
%   p: a vector of single-index Zernike coefficients,
%   waveFront_deltas: wavefronts corresponding to phase diversity images
%   (N)
%   Rc: Rc matrix
%   penalChioce: penalty term to image, 
%       1 for L2 norm of image and 2 for L2 norm of gradient
%   gamma: parameter for the regularization to image
%   alpha: parameter for the regularization to phase
%   itLimit: maximum iteration number
%   stopChoice: iteration stopping criterion
%       0: no stopping criterion;
%       1: based on RMS of wavefront;
%       2: based on loss function;
%   tolValue = 0.01; % iteration stopping criterion: tolerance value
%   GPUflag: GPU options, 0: CPU; 1: GPU;
%   nFlag: Zernike normalization option

% By: Min Guo
% Jan 29, 2020
% Modifications: July 1, 2020
%   add option for L2 norm of gradient: penalChioce = 1 or 2; 
%       1 for L2 norm to image and 2 for L2 norm to gradient
%   modify input pupil parameters(r, theta) from 1D vectors to 2D matrix 
% Modifications: July 18, 2020
%   recover the option of alpha
%   add termination criterion for the iteration


r = r0(idx); % convert 2D matrix to vector with selected elements
theta = theta0(idx); % convert 2D matrix to vector with selected elements

% alpha = 0;
[Sx, Sy, imgNum] = size(imgDs); % image size and numbers
opNum = length(p);
idxNum = length(r);
J = zeros(itLimit+1, 1, 'single'); % cost function values (initial + all iterations)
cInter = zeros(itLimit, opNum, 'single');% Zernike coefficients (all iterations)
Jrms = zeros(itLimit, 1, 'single');
% GPU Data transfering
if(GPUflag==1)
    imgInter = zeros(Sx,Sy,itLimit, 'single', 'gpuArray');
    pupilMask = zeros(Sx, Sy, 'single', 'gpuArray');
    waveFront = zeros(Sx,Sy, 'single', 'gpuArray');
    Hks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    hks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    Sks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    cEstimate = gpuArray(c0);
    cDelta = zeros(size(c0),'single', 'gpuArray');
    zernPolynomial2D = zeros(Sx, Sy, 'single', 'gpuArray');
    
    g = zeros(opNum,1, 'single', 'gpuArray');
    DQks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    HGN_phi_ns = zeros(Sx,Sy,opNum, 'single', 'gpuArray');
    Hmatrix = zeros(opNum,opNum, 'single', 'gpuArray');
    J = gpuArray(J);
    Jrms = gpuArray(Jrms);
    cInter = gpuArray(cInter);
else
    imgInter = zeros(Sx,Sy,itLimit, 'single');
    pupilMask = zeros(Sx, Sy, 'single');
    waveFront = zeros(Sx,Sy, 'single');
    
    Hks = zeros(Sx, Sy, imgNum, 'single');
    hks = zeros(Sx, Sy, imgNum, 'single');
    Sks = zeros(Sx, Sy, imgNum, 'single');
    
    cEstimate = c0;
    cDelta = zeros(size(c0),'single');
    zernPolynomial2D = zeros(Sx, Sy, 'single');
    g = zeros(opNum,1, 'single');
    DQks = zeros(Sx, Sy, imgNum, 'single');
    HGN_phi_ns = zeros(Sx,Sy,opNum, 'single');
    Hmatrix = zeros(opNum,opNum, 'single');
    
end
pupilMask(idx) = 1;

% calculate penalty term
switch penalChioce
    case 1
        penalTerm = gamma;
    case 2
        penalTerm = gamma* r0.*r0; % gamma*u^2
    otherwise
        error('zernretrieve_loop: wrong penalty choice');
end

% for calculating cost function
imgDsSq = sum(abs(imgDs).^2,3); 

for i = 1:itLimit+1
    % ================== %
    % Gradient: g
    % ================== %
    c = cEstimate;
    waveFront(idx) = create_wavefront(p, c, r, theta, nFlag);
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
    
    % calculate Q
    Q = sum(abs(Sks).^2,3) + penalTerm;
    
    % % calculate cost function value: last iteration
    phaseTerm = alpha*c'*Rc*c;
    numeTerm = sum(conj(imgDs).*Sks,3);
    tTerm = imgDsSq - (abs(numeTerm).^2)./Q;
    J(i) = sum(tTerm(:)) + phaseTerm;
    
    % % % % check if termanate iteration
    %       0: no stopping criterion;
    %       1: based on RMS of wavefront;
    %       2: based on loss function;
    switch(stopChoice)
        case 1 % based on RMS of wavefront;
            waveFrontTemp = zernPolynomials*c;
            Jrms(i) = rms(waveFrontTemp);
            waveFrontTemp = zernPolynomials*cDelta;
            rmsDelta = rms(waveFrontTemp);
            if(i>=3)
                if(rmsDelta >= Jrms(i)||rmsDelta >= Jrms(i-1))
                    itTotal = i-2; % diverges, use second last estimate
                    break;
                end
                if(rmsDelta/Jrms(i)<tolValue)
                    itTotal = i-1; % converges, use last estimate
                    break;
                end
            end
        case 2 % based on RMS of wavefront;
            Jlast = J(i);
            if(i>=3)
%               % compare to last two iterations to avoid fluctuations
                if((J(i-2)<Jlast)&&(J(i-1)<Jlast)) % diverges, use second last estimate
                    itTotal = i-2;
                    break;   
                end
            end
            if(i>=2)
                Jdelta = abs(J(i-1) - Jlast)/abs(J(1)-Jlast);
                if(Jdelta<tolValue)
                    itTotal = i-1; % converges, use last estimate
                    break;
                end
            end
        otherwise % no stopping criterion;
    end   
    % if reach maximum iteration
    if(i==itLimit+1)
        itTotal = itLimit;
        break;
    end
    
    % % % % continue iterating
    % calculate F
    numeTerm = sum(conj(Sks).*imgDs,3);
    F = numeTerm./Q;
    
    % calculate Vk and g[phi]
    g_phi = zeros(idxNum,1, 'single');
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
    cDelta = - Hmatrix\g;
    cEstimate = c + cDelta;
    cInter(i, :) = cEstimate;
    imgEstimate = real(ifft2(F));
    imgEstimate(imgEstimate<0) = 0;
    imgInter(:,:,i) = imgEstimate;
             
end
cEstimate = cInter(itTotal, :);
imgEstimate = imgInter(:,:,itTotal);
if (GPUflag == 1)
    J = gather(J);
    cEstimate = gather(cEstimate);
    imgEstimate = gather(imgEstimate);
    cInter = gather(cInter);
end


