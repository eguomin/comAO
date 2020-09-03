function [img3D, PSF3D] = gen_simu_3Dimage(img0, p, coeffs, pixelSize, ...
    lambda, NA, zStepSize, RI, nType)
% simulate 3D image based on ground truth image and zernike coefficients
%
% Output
%   img3D: 3D output images
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

% By: Min Guo
% Apr 10, 2020
[Sx, Sy, Sz] = size(img0);
if(Sx~=Sy)
   error('gen_simu_images: the x size of input images should be same with the y size');
end 

% genearate PSF
PSF3D = coeffs2PSF(p,coeffs, Sx, pixelSize, lambda, NA, Sz,zStepSize, RI);

OTF = genOTF_MATLAB(PSF3D);

imgBlur = ConvFFT3_S(img0, OTF);

switch nType
    case 'none'
        img3D = imgBlur;
    case 'gaussian'
        img3D = zeros(Sx,Sy,Sz, 'single');
        for i = 1:Sz
            img3D(:,:,i) = addgaussiannoise(img3D(:,:,i),0.05);
        end
    case 'poisson'
        img3D = zeros(Sx,Sy,Sz, 'single');
        for i = 1:Sz
            img3D(:,:,i) = addpoissonnoise(img3D(:,:,i));
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
OTF = fftn(ifftshift(PSF));
end

function outVol = ConvFFT3_S(inVol,OTF)
    outVol = real(ifftn(fftn(inVol).*OTF));  
end

