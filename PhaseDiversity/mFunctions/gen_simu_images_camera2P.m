function [imgs, PSFs, waveFronts, PSF_Exc2P_norms] = gen_simu_images_camera2P(img0, p, coeffs, ...
    pixelSize, lambdaExc, lambdaDet, NA, nType, SNR)
% simulate 2P phase diversity images based on ground truth image and phases
%
% Output
%   imgs: phase diversity images, third dimension: N, the number of images
%   waveFronts: wavefront images (in um unit) corresponding to input coeffs
%           The first image is the base or offset coefficients
% Input
%   img0: input ground truth image
%   p: a vector of single indexes(OSA/ANSI convention) for Zernike components,
%   coeffs : a matrix (unknown and known phases) of Zernike coefficients; the
%   first dimension corresponding to p;
%   the second dimension should be the image number N; first row is the
%   base or offset coefficients
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   nType: noise type
%       1)'none': noise free
%       2)'gaussian': Gaussian noise
%       3)'poisson': Poisson noise
%       3)'both': Gaussion and Poisson noise
%   SNR: SNR
%   rotAng: rotation angle between wavefront (sensor) and diveristy 
%       images, (unit: degree) 

% By: Min Guo
% Jan 30, 2020
% Modified: June 19, 2020

% Modified: Aug 29, 2020
% define SNR based on aberrated image (previously based on ground truth)
% Choice for SNR
SNRchoice = 2; % 1: define SNR based on the signal of ground truth
         % 2: define SNR based on the signal of aberrated image (1st image)
% if not define SNR (to be compatible to old version);
if(nargin == 8)
    SNR = 10;
end

[Sx, Sy] = size(img0);
imgNum = size(coeffs,1);
normFlag = 0;
if(Sx~=Sy)
   error('gen_simu_images: the x size of input images should be same with the y size');
end 

if(~strcmp(nType, 'none'))
    % desired signal level
    G = 10; % background white noise (Gaussian noise)
    S = addnoise_snr2signal(SNR, nType, G);
    
    % singal defined based on ground truth image
    if(SNRchoice==1)
        sType = 3; %defined as the intensity average of the pixels > threshold.
        threPerc = 0.1; % threshold value
        S0 = addnoise_getSignal(img0,sType, threPerc);
    end
end

imgFFT = fft2(img0);
imgs = zeros(Sx, Sy, imgNum,'single');
PSFs = zeros(Sx, Sy, imgNum,'single');
waveFronts = zeros(Sx, Sy, imgNum,'single');
PSF_Exc2P_norms = zeros(1,imgNum);
for i = 1:imgNum
    if(i==1)
        iCoeffs = coeffs(1,:);
    else
        iCoeffs = coeffs(1,:) + coeffs(i,:);
    end
    [PSF, ~, ~, ~, PSF_Exc2P_norm]= genPSF_camera2P_2D(p,iCoeffs, Sx, ...
    pixelSize, lambdaExc,lambdaDet, NA, normFlag);

    
    OTF = fft2(ifftshift(PSF));
    imgBlurFFT = imgFFT.* OTF;
    imgBlur = real(ifft2(imgBlurFFT));
    if(i==1)
        % singal defined based on aberrated image (1st image)
        if(SNRchoice==2)
            sType = 3; %defined as the intensity average of the pixels > threshold.
            threPerc = 0.1; % threshold value
            S0 = addnoise_getSignal(imgBlur,sType, threPerc);
        end        
    end
    if(strcmp(nType, 'none'))
        img = imgBlur; % noise free
    else % scale all images by same value and add noise
        % image scaling ratio
        scaleR = S/S0;
        img = addnoiseScale(imgBlur, nType, scaleR, G);
    end
    
    imgs(:,:,i) = img;
    PSFs(:,:,i) = PSF;
    PSF_Exc2P_norms(1,i) = PSF_Exc2P_norm;
end

for i = 1: imgNum
    [waveFront,~, ~, ~] = coeffs2wavefront(p,coeffs(i,:),Sx,...
        pixelSize, lambdaDet, NA);
    waveFronts(:,:,i) = waveFront;
end