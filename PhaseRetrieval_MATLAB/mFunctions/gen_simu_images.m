function [imgs, PSFs, waveFronts] = gen_simu_images(img0, p, coeffs, pixelSize, lambda, NA, nType)
% simulate phase diversity images based on ground truth image and phases

% Output
%   imgs: phase diversity images, third dimension: N, the number of images
%   waveFronts: wavefront images (in um unit) corresponding to input coeffs
%           The first image is the base or offset coefficients
% Input
%   img0: input ground truth image
%   p: a vector of single indexes(Fringe convention) for Zernike components,
%   coeffs : a matrix (unknown and known phases) of Zernike coefficients; the
%   first dimension corresponding to p (should be normalized to phase unit: pi);
%   the second dimension should be the image number N; first row is the
%   base or offset coefficients
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   nType: noise type
%       1)'none': noise free
%       2)'gaussian': Gaussian noise
%       3)'poisson': Poisson noise

% By: Min Guo
% Jan 30, 2020
[Sx, Sy] = size(img0);
imgNum = size(coeffs,1);
if(Sx~=Sy)
   error('gen_simu_images: the x size of input images should be same with the y size');
end 

% Zernike coefficients: convert lengh unit(um) to phase unit(pi)
length2phase = 2*pi/lambda;
coeffs = length2phase * coeffs;

% Define the pupil coordinates (Polar coordinate system) 
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
r0 = r(idx);
theta0 = theta(idx);

phi = zeros(Sx, Sy);
pupilMask = zeros(Sx, Sy);
pupilMask(idx) = 1;
imgFFT = fft2(img0);
imgs = zeros(Sx, Sy, imgNum);
PSFs = zeros(Sx, Sy, imgNum);
waveFronts = zeros(Sx, Sy, imgNum);
for i = 1:imgNum
    if(i==1)
        c_phi = coeffs(1,:);
    else
        c_phi = coeffs(1,:) + coeffs(i,:);
    end
    phi(idx) = create_wavefront(p,c_phi,r0,theta0); % in phase unit: pi
    pupilFun = pupilMask.*exp(1i*phi);
    prf = fftshift(ifft2(ifftshift(pupilFun)));
    PSF = abs(prf).^2;
    OTF = fft2(ifftshift(PSF));
    imgBlurFFT = imgFFT.* OTF;
    imgBlur = ifft2(imgBlurFFT);
    switch nType
        case 'none'
            img = imgBlur;
        case 'gaussian'
            img = addgaussiannoise(imgBlur,0.05);
        case 'poisson'
            img = addpoissonnoise(imgBlur);
        otherwise
            error('gen_simu_images: wrong noise type');
    end
    imgs(:,:,i) = img;
    PSFs(:,:,i) = PSF;
    waveFronts(:,:,i) = phi;
end

% convert phase unit(pi) to lengh unit(um)
waveFronts = waveFronts/length2phase;
% remove offset
if(imgNum > 1)
    for i = 2: imgNum
        waveFronts(:,:,i) = waveFronts(:,:,i) - waveFronts(:,:,1);
    end
end