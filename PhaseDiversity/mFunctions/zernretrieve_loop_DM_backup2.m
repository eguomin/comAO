function [cEstimate, bEstimate, imgEstimate, J, itTotal, cInter, bInter] = zernretrieve_loop_DM(imgDs, ...
    r0, theta0, idx, p, deltaVol, c0, b0, zernPolynomials, Rc, ...
    penalChioce, gamma, alpha, itLimit, stopChoice, tolValue, GPUflag, nFlag)
% retrieve the zernike coefficients of of wavefront slope of an DM actuator 
% based on phase diversity images;
% update formula(recurrence): 
%                       [c_i+1; b_i+1] = [c_i; b_i] + [c_delta; b_delta]
%                       [c_delta; b_delta] = -H^(-1)g
%                       g = [g_c; g_b]
%                       H = [H_cc, g_cb]
%                           [g_bc, H_bb]
% Output
%   cEstimate: retrieved Zernike coefficients(um), initial wavefront
%   bEstimate: retrieved Zernike coefficients(um), delta slope
%   imgEstimate: retrieved image
%   J: cost function values at each iteration
%   itTotal: actual iteration number
%   cInter: Zernike coefficients at each iteration, initial wavefront
%   bInter: Zernike coefficients at each iteration, delta slope
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
% Mar 24, 2022


r = r0(idx); % convert 2D matrix to vector with selected elements
theta = theta0(idx); % convert 2D matrix to vector with selected elements

% alpha = 0;
[Sx, Sy, imgNum] = size(imgDs); % image size and numbers
kStepNum = (imgNum-1)/2;
% kSteps = [0 -kStepNum:-1 1:kStepNum]; % organized to match previous images format
kSteps = -kStepNum:kStepNum; % organized to match previous images format
opNum = length(p);
idxNum = length(r);
J = zeros(itLimit+1, 1, 'single'); % cost function values (initial + all iterations)
cInter = zeros(itLimit, opNum, 'single');% Zernike coefficients (all iterations)
bInter = zeros(itLimit, opNum, 'single');% Zernike coefficients (all iterations)
Jrms = zeros(itLimit, 1, 'single');
% GPU Data transfering
if(GPUflag==1)
    imgInter = zeros(Sx,Sy,itLimit, 'single', 'gpuArray');
    pupilMask = zeros(Sx, Sy, 'single', 'gpuArray');
    waveFront = zeros(Sx,Sy, 'single', 'gpuArray');
    waveFrontb = zeros(Sx,Sy, 'single', 'gpuArray');
    Hks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    hks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    Sks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    cEstimate = gpuArray(c0);
    bEstimate = gpuArray(b0);
    cDelta = zeros(size(c0),'single', 'gpuArray');
    bDelta = zeros(size(b0),'single', 'gpuArray');
    zernPolynomial2D = zeros(Sx, Sy, 'single', 'gpuArray');
    
    % g = zeros(2*opNum,1, 'single', 'gpuArray');
    g_c = zeros(opNum,1, 'single', 'gpuArray');
    g_b = zeros(opNum,1, 'single', 'gpuArray');
    DQks = zeros(Sx, Sy, imgNum, 'single', 'gpuArray');
    HGN_cc_phi_ns = zeros(Sx,Sy,opNum, 'single', 'gpuArray');
    HGN_cb_phi_ns = zeros(Sx,Sy,opNum, 'single', 'gpuArray');
    HGN_bb_phi_ns = zeros(Sx,Sy,opNum, 'single', 'gpuArray');
    Hmatrix = zeros(opNum,opNum, 'single', 'gpuArray');
    J = gpuArray(J);
    Jrms = gpuArray(Jrms);
    cInter = gpuArray(cInter);
    bInter = gpuArray(bInter);
else
    imgInter = zeros(Sx,Sy,itLimit, 'single');
    pupilMask = zeros(Sx, Sy, 'single');
    waveFront = zeros(Sx,Sy, 'single');
    waveFrontb = zeros(Sx,Sy, 'single');
    
    Hks = zeros(Sx, Sy, imgNum, 'single');
    hks = zeros(Sx, Sy, imgNum, 'single');
    Sks = zeros(Sx, Sy, imgNum, 'single');
    
    cEstimate = c0;
    bEstimate = b0;
    cDelta = zeros(size(c0),'single');
    bDelta = zeros(size(b0),'single');
    zernPolynomial2D = zeros(Sx, Sy, 'single');
    % g = zeros(2*opNum,1, 'single');
    g_c = zeros(opNum,1, 'single');
    g_b = zeros(opNum,1, 'single');
    DQks = zeros(Sx, Sy, imgNum, 'single');
    HGN_cc_phi_ns = zeros(Sx,Sy,opNum, 'single');
    HGN_cb_phi_ns = zeros(Sx,Sy,opNum, 'single');
    HGN_bb_phi_ns = zeros(Sx,Sy,opNum, 'single');
    Hmatrix = zeros(2*opNum,2*opNum, 'single');
    Hmatrix_cc = zeros(opNum,opNum, 'single');
    Hmatrix_cb = zeros(opNum,opNum, 'single');
    Hmatrix_bb = zeros(opNum,opNum, 'single');
    
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
    b = bEstimate;
    waveFront(idx) = create_wavefront(p, c, r, theta, nFlag);
    waveFrontb(idx) = create_wavefront(p, b, r, theta, nFlag);
    for k = 1:imgNum
        phi = waveFront+ kSteps(k)* deltaVol *waveFrontb; %phase
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
    phaseTerm = alpha*c'*Rc*c + alpha*b'*Rc*b;
    numeTerm = sum(conj(imgDs).*Sks,3);
    tTerm = imgDsSq - (abs(numeTerm).^2)./Q;
    J(i) = sum(tTerm(:)) + phaseTerm;
    
    % % % % check if termanate iteration
    %       0: no stopping criterion;
    %       1: based on RMS of wavefront;
    %       2: based on loss function;
    switch(stopChoice)
        case 1 % based on RMS of wavefront;
            waveFrontTemp = zernPolynomials*b;
            Jrms(i) = rms(waveFrontTemp);
            waveFrontTemp = zernPolynomials*bDelta;
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
    if(GPUflag==1)
        g_c_phi = zeros(idxNum,1, 'single', 'gpuArray');
        g_b_phi = zeros(idxNum,1, 'single', 'gpuArray');
    else
        g_c_phi = zeros(idxNum,1, 'single');
        g_b_phi = zeros(idxNum,1, 'single');
    end
    for k = 1:imgNum
        Sk = squeeze(Sks(:,:,k));
        Dk = squeeze(imgDs(:,:,k));
        Vk = conj(F).*Dk - abs(F).^2.*Sk;
        Vk_iFT = ifft2(Vk);
        hk = hks(:,:,k);
        Hk = Hks(:,:,k);
        temp1 = hk.*real(Vk_iFT);
        % temp1 = real(conj(hk)).*Vk_iFT;
        temp2 = fft2(temp1);
        temp3 = imag(conj(Hk).*temp2);
        g_c_phi = g_c_phi - 2*temp3(idx);
        g_b_phi = g_b_phi - 2*kSteps(k)*deltaVol*temp3(idx);
    end
    
    % calculate gradient: g
    for m = 1:opNum
        zernPolynomial = zernPolynomials(:,m);
        g_c(m) = dot(g_c_phi,zernPolynomial);
        g_b(m) = dot(g_b_phi,zernPolynomial);
    end
    g_c = g_c + 0.5 * alpha * Rc * c;
    g_b = g_c + 0.5 * alpha * Rc * c;
    
    g = [g_c; g_b];
    
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
        HGN_cc_phi_n = 0;
        HGN_cb_phi_n = 0;
        % HGN_bc_phi_n = 0;
        HGN_bb_phi_n = 0;
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
                % temp1 = ifft2(conj(DQk).*conj(Ujk));
                temp2 = fft2(hj.*temp1);
                temp3 = conj(Hj).*temp2;
            
                temp4 = ifft2(conj(DQj).*Ujk);
                temp5 = fft2(hk.*temp4);
                temp6 = conj(Hk).*temp5;
                
                kd = kSteps(k)*deltaVol;
                HGN_cc_phi_n = 4* imag(temp3-temp6) + HGN_cc_phi_n; 
                HGN_cb_phi_n = 4* kd*imag(temp3-temp6) + HGN_cb_phi_n;
                % HGN_bc_phi_n = HGN_cb_phi_n;
                HGN_bb_phi_n = 4* kd*kd*imag(temp3-temp6) + HGN_bb_phi_n;
            end
        end
        HGN_cc_phi_ns(:,:,n) = HGN_cc_phi_n;
        HGN_cb_phi_ns(:,:,n) = HGN_cb_phi_n;
        HGN_bb_phi_ns(:,:,n) = HGN_bb_phi_n;
    end
        
    % calculate Hessian Matrix elements
    for m = 1:opNum
        zernPolynomial = zernPolynomials(:,m); % vector corresponding to idx
        for n = 1:opNum
            HGN_phi_n = squeeze(HGN_cc_phi_ns(:,:,n));
            HGN_phi_n_idx = HGN_phi_n(idx);
            Hmatrix_cc(m,n) = dot(HGN_phi_n_idx, zernPolynomial);
            
            HGN_phi_n = squeeze(HGN_cb_phi_ns(:,:,n));
            HGN_phi_n_idx = HGN_phi_n(idx);
            Hmatrix_cb(m,n) = dot(HGN_phi_n_idx, zernPolynomial);
            
            HGN_phi_n = squeeze(HGN_bb_phi_ns(:,:,n));
            HGN_phi_n_idx = HGN_phi_n(idx);
            Hmatrix_bb(m,n) = dot(HGN_phi_n_idx, zernPolynomial);
        end
    end
    
    Hmatrix(1:opNum,1:opNum) = Hmatrix_cc + alpha*Rc;
    Hmatrix(1:opNum,opNum+1:2*opNum) = Hmatrix_cb + alpha*Rc;
    Hmatrix(opNum+1:2*opNum,1:opNum) = Hmatrix_cb + alpha*Rc;
    Hmatrix(opNum+1:2*opNum,opNum+1:2*opNum) = Hmatrix_bb + alpha*Rc;
    
    combDelta = - Hmatrix\g;
    cDelta = combDelta(1:opNum);
    cEstimate = c + cDelta;
    bDelta = combDelta(opNum+1:2*opNum);
    bEstimate = b + bDelta;
%     c_show = gather(c')
%     cDelta_show = gather(cDelta')
%     cEstimate_show = gather(cEstimate')
%     b_show = gather(b')
%     bDelta_show = gather(bDelta')
%     bEstimate_show = gather(bEstimate')
    cInter(i, :) = cEstimate;
    bInter(i, :) = bEstimate;
    imgEstimate = real(ifft2(F));
    imgEstimate(imgEstimate<0) = 0;
    imgInter(:,:,i) = imgEstimate;
             
end
cEstimate = cInter(itTotal, :);
bEstimate = bInter(itTotal, :);
imgEstimate = imgInter(:,:,itTotal);
if (GPUflag == 1)
    J = gather(J);
    cEstimate = gather(cEstimate);
    bEstimate = gather(bEstimate);
    imgEstimate = gather(imgEstimate);
    cInter = gather(cInter);
    bInter = gather(bInter);
end


