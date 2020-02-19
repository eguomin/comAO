function [coeffsOut, cTime] = recon_zern_simuBatch(pIn,coeffsIn, img0, coeffs_delta, ...
    zernCoeffOrder, iteration, gamma, pixelSize, lambda, NA, nType, flagGPU)

flagSmoothEdge = 0;
if(zernCoeffOrder < max(pIn))
    error('Zernike coefficents order is more than input orders');
end
pEnd = zernCoeffOrder - pIn(1) + 1;
[expNum, zernNum] = size(coeffsIn);
[pdNum, zernNum2] = size(coeffs_delta);
if(zernNum~=zernNum2)
    error('Zernike coefficents number does not match');
end
imgNum = pdNum + 1;
tic;
% upboosting steps,e.g., 4th (15), 5th(21), 6th(28) and 7th(36) orders
if(zernCoeffOrder<=11)
    zernSteps = pIn(zernNum) - pIn(1) + 1;
elseif(zernCoeffOrder<=15)
    zernSteps = [11 pIn(zernNum)] - pIn(1) + 1; 
elseif(zernCoeffOrder<=21)
    zernSteps = [11 15 pIn(zernNum)] - pIn(1) + 1; 
elseif(zernCoeffOrder<=28)
    zernSteps = [11 15 21 pIn(zernNum)] - pIn(1) + 1; 
elseif(zernCoeffOrder<=36)
    zernSteps = [11 15 21 28 pIn(zernNum)] - pIn(1) + 1; 
else
    error('Zernike coefficents order out of range');
end
% pixelSize = 0.096; % um
% lambda = 0.550; % um
% NA = 1.2;
pIn0 = pIn(1:pEnd);
coeffsOut = zeros(expNum, zernNum);
for i = 1:expNum
    cTime1 = toc;
    disp(['Experiment #: ', num2str(i)]);
    coeffsInitial = coeffsIn(i,:);
    coeffs_all = [coeffsInitial;coeffs_delta;];
    % generate phase diversity images;
    disp('...Generating simulated images...');
    [imgs, ~, ~] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, lambda, NA, nType);
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
    [cEstimate, ~, ~] = recon_zern(imgs, pIn0, coeffs_delta0, gamma, iteration, zernSteps, pixelSize, lambda, NA, flagGPU);
    cTime2 = toc;
    coeffsOut(i,1:pEnd) = cEstimate;
    disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
end
cTime = toc;
disp(['Simulation completed!!! Total time cost:', num2str(cTime)]);
