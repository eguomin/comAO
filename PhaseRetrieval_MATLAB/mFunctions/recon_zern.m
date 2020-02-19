function [cEstimate, imgEstimate, rePar] = recon_zern(imgs, p, c_delta, gamma, it, zernSteps, pixelSize, lambda, NA, GPUflag)
% reconstruct the zernike coefficients of aberrated phase 
% based on phase diversity images and upboosting strategy

% Output
%   cEstimate: reconstructed zernike coefficients(um)
%   imgEstimate: estimated image
%   rePar: records of intermidiate parameters
% Input
%   imgs: phase diversity images, third dimension: N, the number of images
%   p: a vector of single indexes(Fringe convention, first coefficient is piston) for Zernike components,
%   c_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (unit: um); if matrix, the second dimension should be N-1
%   gamma: parameter for the regularization
%   it: iteration number
%   zernSteps: coefficients in p to be updated step by step
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   GPUflag: GPU options, 0: CPU; 1: GPU;
% By: Min Guo
% Jan 28, 2020
[Sx, Sy, imgNum] = size(imgs);
deltaNum = size(c_delta,1);
if(Sx~=Sy)
   error('recon_zern: the x size of input images should be same with the y size');
end 
if (imgNum - deltaNum)~= 1
    error('recon_zern:NMlength','deltaNum should be: imgNum -1.')
end
if(length(p)<max(zernSteps(:)))
    error('recon_zern: upSteps exceed the p index');
end

% Zernike coefficients: convert lengh unit(um) to phase unit(pi)
length2phase = 2*pi/lambda;
c_delta = length2phase * c_delta;

% Define the pupil coordinates (Polar coordinate system) 
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
r0 = r(idx);
theta0 = theta(idx);

% pre-calculation
[imgDs, zernPolynomials, Rc, waveFront_deltas] = zernretrieve_pre(imgs, r0, theta0, idx, p, c_delta);

% reconstruction
stepNum = length(zernSteps);
rePar.cEstimate = zeros(length(p),stepNum);
rePar.J = zeros(it,stepNum);
cEstimate = zeros(zernSteps(1),1);
if(GPUflag == 1)
    imgDs = gpuArray(imgDs);
    zernPolynomials = gpuArray(zernPolynomials);
    Rc = gpuArray(Rc);
    waveFront_deltas = gpuArray(waveFront_deltas);
end
for iStep = 1:stepNum
    zernN = zernSteps(iStep);
    p0 = p(1:zernN);
    c0 = zeros(zernN,1);
    zernNLastIt = length(cEstimate);
    c0(1:zernNLastIt) = cEstimate;
    zernPolynomials0 = zernPolynomials(:,1:zernN);
    Rc0 = Rc(1:zernN,1:zernN);
    [cEstimate, imgEstimate, J] = zernretrieve_loop(imgDs, r0, theta0, idx,...
    p0, waveFront_deltas, c0, zernPolynomials0, Rc0, gamma, it, GPUflag);
    rePar.cEstimate(1:zernN,iStep) = cEstimate;
    rePar.J(:,iStep) = J;
end

% convert phase unit(pi) to lengh unit(um)
cEstimate = cEstimate/length2phase;
rePar.cEstimate = rePar.cEstimate/length2phase;