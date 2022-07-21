function [cEstimate, bEstimate, imgEstimate, rePar] = recon_zern_DM(imgs, p, deltaVol, gamma,...
    alpha, itLimit, zernSteps, pixelSize, lambda, NA, GPUflag, penalChioce,rotAng, bInput)
% reconstruct the zernike coefficients of wavefront slope of an DM actuator 
% based on phase diversity images and upboosting strategy
% 
% Output
%   cEstimate: reconstructed zernike coefficients(um)
%   imgEstimate: estimated image
%   rePar: records of intermidiate parameters
% Input
%   imgs: phase diversity images, third dimension: N, the number of images
%   p: a vector of single-index Zernike coefficients,
%   deltaVol: delta voltage on actuator
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
% Mar 24, 2022

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
if(Sx~=Sy)
   error('recon_zern: the x size of input images should be same with the y size');
end 

if(length(p)<max(zernSteps(:)))
    error('recon_zern: upSteps exceed the p index');
end

% Zernike coefficients: convert lengh unit(um) to phase unit(pi)
length2phase = 2*pi/lambda;
c_delta = zeros(imgNum-1, length(p), 'single');
bInput = bInput*length2phase;

% Define the pupil coordinates (Polar coordinate system) 
[r0, theta0, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
theta0 = theta0 + deg2rad(rotAng);

% pre-calculation
[imgDs, zernPolynomials, Rc, ~] = zernretrieve_pre(imgs, ...
    r0, theta0, idx, p, c_delta, nFlag);

% reconstruction
stepNum = length(zernSteps);
rePar.cEstimate = zeros(length(p),stepNum, 'single');
rePar.bEstimate = zeros(length(p),stepNum, 'single');
rePar.J = zeros(itLimit+1,stepNum, 'single');
rePar.itTotal = zeros(1,stepNum, 'single');
rePar.cInter = cell(1,stepNum);
rePar.bInter = cell(1,stepNum);
% cEstimate = zeros(zernSteps(1),1, 'single');
% bEstimate = zeros(zernSteps(1),1, 'single');
cEstimate = 0.1*rand(zernSteps(1),1);
bEstimate = 0.1*rand(zernSteps(1),1);
if(GPUflag == 1)
    imgDs = gpuArray(imgDs);
    zernPolynomials = gpuArray(zernPolynomials);
    Rc = gpuArray(Rc);
end
for iStep = 1:stepNum
    zernN = zernSteps(iStep);
    zernNLastIt = length(cEstimate);
    p0 = p(1:zernN);
%     c0 = zeros(zernN,1, 'single');
    c0 = 0.1*rand(zernN,1);
%     c0(1:zernNLastIt) = cEstimate;
    c0(1:zernN) = bInput(1:zernN);
%     b0 = zeros(zernN,1, 'single');
    b0 = 0.1*rand(zernN,1);
    b0(1:zernNLastIt) = bEstimate;
%     b0(1:zernN) = bInput(1:zernN);
    zernPolynomials0 = zernPolynomials(:,1:zernN);
    Rc0 = Rc(1:zernN,1:zernN);
    [cEstimate, bEstimate, imgEstimate, J, itTotal, cInter, bInter] = zernretrieve_loop_DM(imgDs,...
        r0, theta0, idx, p0, deltaVol, c0, b0, zernPolynomials0, Rc0,...
        penalChioce, gamma, alpha, itLimit, stopChoice, tolValue, GPUflag, nFlag);
    rePar.cEstimate(1:zernN,iStep) = cEstimate;
    rePar.bEstimate(1:zernN,iStep) = bEstimate;
    rePar.J(:,iStep) = J;
    rePar.itTotal(iStep) = itTotal;
    rePar.cInter{iStep} = cInter;
    rePar.bInter{iStep} = bInter;
end

% convert phase unit(pi) to lengh unit(um)
cEstimate = cEstimate/length2phase;
rePar.cEstimate = rePar.cEstimate/length2phase;
bEstimate = bEstimate/length2phase;
rePar.bEstimate = rePar.bEstimate/length2phase;