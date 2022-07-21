function [cEstimate, imgEstimate, rePar] = recon_zern(imgs, p, c_delta, gamma,...
    alpha, itLimit, zernSteps, pixelSize, lambda, NA, GPUflag, penalChioce,rotAng)
% reconstruct the zernike coefficients of aberrated phase 
% based on phase diversity images and upboosting strategy
% 
% Output
%   cEstimate: reconstructed zernike coefficients(um)
%   imgEstimate: estimated image
%   rePar: records of intermidiate parameters
% Input
%   imgs: phase diversity images, third dimension: N, the number of images
%   p: a vector of single-index Zernike coefficients,
%   c_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (unit: um); if matrix, the second dimension should be N-1
%   gamma: parameter for the regularization to image
%   alpha: parameter for the regularization to phase
%   itLimit: maximum iteration number
%   zernSteps: coefficients in p to be updated step by step
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   GPUflag: GPU options, 0: CPU; 1: GPU;
%   penalChioce: penalty term to image, 
%       1 for L2 norm of image and 2 for L2 norm of gradient
%   rotAng: rotation angle between wavefront (sensor) and diveristy 
%       images, (unit: degree) 
      
% By: Min Guo
% Jan 28, 2020
% Modifications: Feb 24, 2020
%   add option on the normalization of Zernike
%   nFlag = 'norm' or 'none'; % normalized zernike bases or not
% Modifications: July 1, 2020
%   add option for L2 norm to gradient: penalChioce = 1 or 2; 
%       1 for L2 norm to image and 2 for L2 norm to gradient
% Modifications: Aug. 8, 2020
%   add option for correcting rotation between wavefront and images
if(nargin==11)
    penalChioce = 1; % default for Tikhonov norm
    rotAng = 0;
elseif (nargin==12)
    rotAng = 0;
end
% % % customize a few settings:
nFlag = 'none'; % Zernike norm option
% iteration stopping criterion
%       0: no stopping criterion;
%       1: based on RMS of wavefront; typically tolValue = 0.01;
%       2: based on loss function; typically tolValue = 0.001;
stopChoice = 0;
tolValue = 0.001; % tolerance value 

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
[r0, theta0, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
theta0 = theta0 + deg2rad(rotAng);

% pre-calculation
[imgDs, zernPolynomials, Rc, waveFront_deltas] = zernretrieve_pre(imgs, ...
    r0, theta0, idx, p, c_delta, nFlag);

% reconstruction
stepNum = length(zernSteps);
rePar.cEstimate = zeros(length(p),stepNum, 'single');
rePar.J = zeros(itLimit+1,stepNum, 'single');
rePar.itTotal = zeros(1,stepNum, 'single');
rePar.cInter = cell(1,stepNum);
cEstimate = zeros(zernSteps(1),1, 'single');
if(GPUflag == 1)
    imgDs = gpuArray(imgDs);
    zernPolynomials = gpuArray(zernPolynomials);
    Rc = gpuArray(Rc);
    waveFront_deltas = gpuArray(waveFront_deltas);
end
for iStep = 1:stepNum
    zernN = zernSteps(iStep);
    p0 = p(1:zernN);
    c0 = zeros(zernN,1, 'single');
    zernNLastIt = length(cEstimate);
    c0(1:zernNLastIt) = cEstimate;
    zernPolynomials0 = zernPolynomials(:,1:zernN);
    Rc0 = Rc(1:zernN,1:zernN);
    [cEstimate, imgEstimate, J, itTotal, cInter] = zernretrieve_loop(imgDs,...
        r0, theta0, idx, p0, waveFront_deltas, c0, zernPolynomials0, Rc0,...
        penalChioce, gamma, alpha, itLimit, stopChoice, tolValue, GPUflag, nFlag);
    rePar.cEstimate(1:zernN,iStep) = cEstimate;
    rePar.J(:,iStep) = J;
    rePar.itTotal(iStep) = itTotal;
    rePar.cInter{iStep} = cInter;
end

% convert phase unit(pi) to lengh unit(um)
cEstimate = cEstimate/length2phase;
rePar.cEstimate = rePar.cEstimate/length2phase;