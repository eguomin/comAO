function [coeffsOut_TN, coeffsOut_TV, cTime] = recon_zern_simuBatch_TV(pIn,coeffsIn, img0, coeffs_delta, ...
    zernCoeffOrder, iteration, gamma, pixelSize, lambda, NA, nType, SNR, flagGPU)
% reconstruct a series of zernike coefficients based on 
% simulated images from a series of diversity phases 
% Output
%   coeffsOut_TN: reconstructed zernike coefficients(um)
%   cTime: time cost
% Input
%   pIn: a vector of single indexes(OSA/ANSI convention) for Zernike components,
%   coeffsIn: input ground truth zernike coefficents, a series of vectors
%   img0: input ground truth image
%   coeffs_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (unit: um); if matrix, the second dimension should be N-1
%   zernCoeffOrder: the number of coefficients index to be simulated and reconstructed
%   iteration: iteration number
%   gamma: parameter for the regularization
%   pixelSize: pixel size of the intensity image, unit: um
%   lambda: wavelength, unit: um
%   NA: numerical aperture of objective
%   nType: noise type
%       1)'none': noise free
%       2)'gaussian': Gaussian noise
%       3)'poisson': Poisson noise
%   SNR: signal to noise ratio
%   flagGPU: GPU options, 0: CPU; 1: GPU;
% By: Min Guo
% Mar 10, 2020
% Modifications: July 2, 2020
%   add option for L2 norm of gradient: penalChioce = 1 or 2; 
%       1 for L2 norm of image and 2 for L2 norm of gradient

flagSmoothEdge = 0;
if(zernCoeffOrder < max(pIn))
    error('Zernike coefficents order is more than input orders');
end
% SNR = 20;
pEnd = zernCoeffOrder - pIn(1) + 1;
[expNum, zernNum] = size(coeffsIn);
[pdNum, zernNum2] = size(coeffs_delta);
if(zernNum~=zernNum2)
    error('Zernike coefficents number does not match');
end
imgNum = pdNum + 1;
tic;
% upboosting steps,e.g., 3th(9), 4th(14), 5th(20), 6th(27) and 7th(35) orders
if(zernCoeffOrder<=9)
    zernSteps = pIn(zernNum) - pIn(1) + 1;
elseif(zernCoeffOrder<=14)
    zernSteps = [9 pIn(zernNum)] - pIn(1) + 1; 
elseif(zernCoeffOrder<=20)
    zernSteps = [9 14 pIn(zernNum)] - pIn(1) + 1; 
elseif(zernCoeffOrder<=27)
    zernSteps = [9 14 20 pIn(zernNum)] - pIn(1) + 1; 
elseif(zernCoeffOrder<=35)
    zernSteps = [9 14 20 27 pIn(zernNum)] - pIn(1) + 1; 
else
    error('Zernike coefficents order out of range');
end
% pixelSize = 0.096; % um
% lambda = 0.550; % um
% NA = 1.2;
pIn0 = pIn(1:pEnd);
coeffsOut_TN = zeros(expNum, zernNum);
coeffsOut_TV = zeros(expNum, zernNum);
for i = 1:expNum
    cTime1 = toc;
    disp(['Experiment #: ', num2str(i)]);
    coeffsInitial = coeffsIn(i,:);
    coeffs_all = [coeffsInitial;coeffs_delta;];
    % generate phase diversity images;
    disp('...Generating simulated images...');
    [imgs, ~, ~] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, lambda, NA, nType, SNR);
    % smooth the boundary of images
    if(flagSmoothEdge==1)
        %sKernel = fspecial('gaussian',[10 10],3);
        napodize = 20;
        for j = 1:imgNum
            img = imgs(:,:,j);
            %img = edgetaper(img,sKernel);
            img = apodize(img, napodize);
            imgs(:,:,j) = img;
        end
    end
    
    disp('...Reconstructing Zernike coefficients...');
    coeffs_delta0 = coeffs_delta(:,1:pEnd);
    penalChioce = 1;
    [cEstimate, ~, ~] = recon_zern(imgs, pIn0, coeffs_delta0, gamma, iteration,...
        zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);
    coeffsOut_TN(i,1:pEnd) = cEstimate;
    penalChioce = 2;
    [cEstimate, ~, ~] = recon_zern(imgs, pIn0, coeffs_delta0, gamma, iteration,...
        zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);
    coeffsOut_TV(i,1:pEnd) = cEstimate;
    cTime2 = toc;
    disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
end
cTime = toc;
disp(['Simulation completed!!! Total time cost:', num2str(cTime)]);
