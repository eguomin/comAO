function coeffsOut = recon_zern_simuBatch(pIn,coeffsIn, img0, coeffs_delta, ...
    gamma,alpha, itLimit,zernSteps, pixelSize, lambda, NA, nType, SNR, flagGPU, penalChioce)
% simulate wavefront reconstruction in batch mode:
% 1) a ground truth image and a series of ground truth Zernike coefficients
% 2) sreconstruct a series of zernike coefficients based on 
% simulated images from a series of diversity phases 
% Output
%   coeffsOut: reconstructed zernike coefficients(um)
%   cTime: time cost
% Input
%   pIn: a vector of single indexes(OSA/ANSI convention) for Zernike components,
%   coeffsIn: ground truth zernike coeffs, 2D matrix(a series of vectors)
%   img0: ground truth image
%   coeffs_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to pIn (unit: um); if matrix, the second dimension should be N-1
%   gamma: parameter for the regularization to image
%   alpha: parameter for the regularization to phase
%   itLimit: maximum iteration number
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   nType: noise type
%       1)'none': noise free
%       2)'gaussian': Gaussian noise
%       3)'poisson': Poisson noise
%       3)'both': Gaussion and Poisson noise
%   SNR: SNR
%   flagGPU: GPU options, 0: CPU; 1: GPU;
%   penalChioce: penalty term to image, 
%       1 for L2 norm of image and 2 for L2 norm of gradient

% By: Min Guo
% Mar 10, 2020
% Modifications: add penalty choice and change zernSteps
% Aug 15, 2020

flagSmoothEdge = 0;
[expNum, zernNum] = size(coeffsIn);
[pdNum, zernNum2] = size(coeffs_delta);
if(zernNum~=zernNum2)
    error('Zernike coefficents number does not match');
end
imgNum = pdNum + 1;
tStart = tic;
disp('Batch simulation starting...');

coeffsOut = zeros(expNum, zernNum);
for i = 1:expNum
    cTime1 = toc(tStart);
    disp(['Experiment #: ', num2str(i)]);
    coeffsInitial = coeffsIn(i,:);
    coeffs_all = [coeffsInitial;coeffs_delta;];
    % generate phase diversity images;
    disp('...Generating simulated images...');
    [imgs, ~, ~] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, lambda, NA, nType, SNR);
    % smooth the boundary of images
    if(flagSmoothEdge==1)
        sKernel = fspecial('gaussian',[10 10],3);
        for j = 1:imgNum
            img = imgs(:,:,j);
            img = edgetaper(img,sKernel);
            imgs(:,:,j) = img;
        end
    end
    cTime2 = toc(tStart);
    disp(['... ... time cost: ', num2str(cTime2-cTime1), ' s']);
    
    disp('...Reconstructing Zernike coefficients...');
    [cEstimate, ~, ~] = recon_zern(imgs, pIn, coeffs_delta, gamma, alpha, ...
    itLimit,zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);
    cTime3 = toc(tStart);
    coeffsOut(i,:) = cEstimate;
    disp(['... ... time cost: ', num2str(cTime3-cTime2), ' s']);
    disp(['Current experiment time cost: ', num2str(cTime3-cTime1), ' s']);
end
cTime = toc(tStart);
disp(['Batch simulation completed!!! Time cost:', num2str(cTime), ' s']);
