function [img3D, PSF, PSF0] = gen_simu_3Dimage_fromReal(img0, p, coeffs, pixelSize, ...
    lambda, NA, zStepSize, RI, nType, psfChoice, LsFWHMz)
% simulate 3D aberrated wide-field or light-sheet image based on real image
% (aberration free) and zernike coefficients
%
% Output
%   imgOut: 3D output images
%   PSF: aberrated PSF corresponding to input coeffs
%   PSF0: aberration free PSF
% Input
%   img0: input real image (experimental image without aberration)
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
%   psfChoice: PSF type, 0: wide-field PSF; 1: light-sheet PSF; 2: confocal
%   PSF or two-photon PSF (normalized)
%   LsFWHMz: light-sheet thickness, unit: um

% By: Min Guo
% Sep. 28, 2020
[Sx0, Sy0, Sz0] = size(img0);

% psfChoice = 1; % 0: wide-field PSF; 1: light-sheet PSF; 2: confocal PSF
% LsFWHMz = 2.5; % light-sheet thickness, unit: um
% pixelSizez = 1;
flagPad = 1;
padx = 0;
pady = 0;
padz = 128;

if(flagPad==1)
    Sx = Sx0 + padx;
    Sy = Sy0 + pady;
    Sz = Sz0 + padz;
    img = alignsize3d(img0, [Sx, Sy, Sz]);
else
    Sx = Sx0;
    Sy = Sy0;
    Sz = Sz0;
    img = img0;
end

% genearate PSFs
PSFx = 128;
% PSFy = FPSFx;
PSFz = 128;
coeffs0 = zeros(size(coeffs));
PSF0 = coeffs2PSF(p,coeffs0, PSFx, pixelSize, lambda, NA, PSFz,zStepSize, RI);
if(psfChoice==1)
    PSF0 = genPSF_wf2ls(PSF0, LsFWHMz, zStepSize);
end
if(psfChoice==2)
    PSF0 = PSF0.^2;
end
OTF0 = genOTF_MATLAB(PSF0, [Sx, Sy, Sz]);
PSF = coeffs2PSF(p,coeffs, PSFx, pixelSize, lambda, NA, PSFz,zStepSize, RI);
if(psfChoice==1)
    PSF = genPSF_wf2ls(PSF, LsFWHMz, zStepSize);
end
if(psfChoice==2)
    PSF = PSF.^2;
end
OTF = genOTF_MATLAB(PSF, [Sx, Sy, Sz]);

% blur image
thresValue = 0.01;
OTF_temp = abs(OTF0);
OTF0(OTF_temp<thresValue) = 1;

imgFFT = fftn(img);
imgFFT_blurred = imgFFT.*OTF./OTF0;
imgBlur = real(ifftn(imgFFT_blurred));

if(flagPad==1)
    imgBlur = alignsize3d(imgBlur, [Sx0, Sy0, Sz0]);
end
imgBlur = max(imgBlur,0);

switch nType
    case 'none'
        img3D = imgBlur;
    case 'gaussian'
        img3D = zeros(Sx0,Sy0,Sz0, 'single');
        for i = 1:Sz0
            img3D(:,:,i) = addgaussiannoise(imgBlur(:,:,i),0.05);
        end
    case 'poisson'
        img3D = zeros(Sx0,Sy0,Sz0, 'single');
        for i = 1:Sz0
            img3D(:,:,i) = addpoissonnoise(imgBlur(:,:,i));
        end
    otherwise
        error('gen_simu_images: wrong noise type');
end

end

% functions for convovlution in Fourier domain
function OTF = genOTF_MATLAB(PSF, imgSize)
% calculate OTF for matlab deconvolution
% % output
% OTF:
% % input 
% PSF: input images
% imSize: image size [Sx,Sy,Sz]
if(nargin==2)
    PSF = single(alignsize3d(PSF, imgSize));
end
PSF = PSF/sum(PSF(:));
OTF = fftn(ifftshift(PSF));
end
