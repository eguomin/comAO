% Explore graduated optimization scheme based on simulated phase diversity images

% By: Min Guo
% Last update: Sep. 1, 2020

% clear all;
% close all;
% % % default settings
tStart = tic;
flagGPU = 0;
if(flagGPU)
    devNum = 1;
    gpuDevice(devNum);
end
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
nType = 'both';
SNR = 20;

% % % Unknown aberrations: to be estimated 
flagSetCoeffs = 1; % Set input Zernike coefficients:
if(flagSetCoeffs==1)
    coeffsFile = '..\..\..\computAO\Data\abeAndNoise\aberrationRondom2_0p15.mat';
    varNames = cell(1,3);
    varNames{1} = 'pIn';
    varNames{2} = 'aValue';
    varNames{3} = 'coeffs';
    importfile(coeffsFile,varNames);
    zernNum = length(pIn);
    expNum = size(coeffs,1);
elseif(flagSetCoeffs==2)
    expNum = 50;
    zernNum = length(pIn);
    coeffs = zeros(expNum,zernNum);
    aValue = 0.15;
    zernType = 'random2'; % 'defocus','astig', 'coma','trefoil','sphe', 'random1'
    for i = 1:expNum      
        % % % higher weight for lower order indices
        zernType = 'random2';
        c = gen_zern_coeffs(pIn,aValue,zernType);
        coeffs(i,:) = c;
    end
end

% % % diversity images:
imgNum = 2; % 2 or 3
abType1 = 'defocus'; % 'astig' or 'defocus'
abValue1 = 1; % um
flagSmoothEdge = 0;

% % % reconstruction settings
itLimit = 10; % max iteration number
gamma = 1e-10; % % gamma range: 1e-14 ~ 1e-6
alpha = 0;
penalChioce = 1;
zernSteps1 = pIn(zernNum) - pIn(1) + 1;
% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
zernSteps2 = [14 20 27 pIn(zernNum)] - pIn(1) + 1; 

% % % files and folders
fileFolderIn = '..\DataForTest\Simu\';
fileNameIn = 'imSample1';
fileImgSample = [fileFolderIn, fileNameIn, '.tif'];
img0 = single(ReadTifStack(fileImgSample));
[Sx, ~] = size(img0);
fileFolderOut = ['..\..\..\computAO\Data\bootstrapping_batch\',...
    fileNameIn, '_rd2_0p15_defocus_1\'];
fileNameOut = [nType, '_SNR_', num2str(SNR)];
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

% Known aberration: for diverisity image 1
coeffs_delta = gen_zern_coeffs(pIn,abValue1,abType1);

% % % % **** perform simulations in batch ****
% % % noise free case
coeffsEst1 = zeros(expNum, zernNum);
coeffsEst2 = zeros(expNum, zernNum);
for i = 1:expNum
    cTime1 = toc(tStart);
    disp(['Experiment #: ', num2str(i)]);
    coeffsInitial = coeffs(i,:);
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
    
    disp('...direct reconstruction...');
    [cEstimate1, ~, ~] = recon_zern(imgs, pIn, coeffs_delta, gamma, alpha, ...
    itLimit,zernSteps1, pixelSize, lambda, NA, flagGPU, penalChioce);
    coeffsEst1(i,:) = cEstimate1;
    cTime3 = toc(tStart);
    disp(['... ... time cost: ', num2str(cTime3-cTime2), ' s']);
    
    disp('...bootstrapping reconstruction...');
    [cEstimate2, ~, ~] = recon_zern(imgs, pIn, coeffs_delta, gamma, alpha, ...
    itLimit, zernSteps2, pixelSize, lambda, NA, flagGPU, penalChioce);
    coeffsEst2(i,:) = cEstimate2;
    cTime4 = toc(tStart);
    disp(['... ... time cost: ', num2str(cTime4-cTime3), ' s']);
    disp(['Current experiment time cost: ', num2str(cTime4-cTime1), ' s']);
end
% % % statistics
disp('Performing statistics...');
staWaveFrontIn = zeros(expNum,2);
staWaveFrontEst1 = zeros(expNum,2);
staWaveFrontEst2 = zeros(expNum,2);
staWaveFrontEstError1 = zeros(expNum,2);
staWaveFrontEstError2 = zeros(expNum,2);
for i = 1:expNum
    c = coeffs(i,:);
    c1 = coeffsEst1(i,:);
    c2 = coeffsEst2(i,:);
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,c,Sx,...
        pixelSize, lambda, NA, 0);
    staWaveFrontIn(i,1) = staPara.rmsLength;
    staWaveFrontIn(i,2) = staPara.pvLength;
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,c1,Sx,...
        pixelSize, lambda, NA, 0);
    staWaveFrontEst1(i,1) = staPara.rmsLength;
    staWaveFrontEst1(i,2) = staPara.pvLength;
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,c1-c,Sx,...
        pixelSize, lambda, NA, 0);
    staWaveFrontEstError1(i,1) = staPara.rmsLength;
    staWaveFrontEstError1(i,2) = staPara.pvLength;
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,c2,Sx,...
        pixelSize, lambda, NA, 0);
    staWaveFrontEst2(i,1) = staPara.rmsLength;
    staWaveFrontEst2(i,2) = staPara.pvLength;
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,c2-c,Sx,...
        pixelSize, lambda, NA, 0);
    staWaveFrontEstError2(i,1) = staPara.rmsLength;
    staWaveFrontEstError2(i,2) = staPara.pvLength;
end
staWaveFrontIn_Ave = mean(staWaveFrontIn, 1);
staWaveFrontIn_SD = sqrt(var(staWaveFrontIn, 1));
staWaveFrontEst1_Ave = mean(staWaveFrontEst1, 1);
staWaveFrontEst1_SD = sqrt(var(staWaveFrontEst1, 1));
staWaveFrontEstError1_Ave = mean(staWaveFrontEstError1, 1);
staWaveFrontEstError1_SD = sqrt(var(staWaveFrontEstError1, 1));
staWaveFrontEst2_Ave = mean(staWaveFrontEst2, 1);
staWaveFrontEst2_SD = sqrt(var(staWaveFrontEst2, 1));
staWaveFrontEstError2_Ave = mean(staWaveFrontEstError2, 1);
staWaveFrontEstError2_SD = sqrt(var(staWaveFrontEstError2, 1));
cTime = toc(tStart);
save([fileFolderOut, fileNameOut,'.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime), ' s']);

% % % Check 1 reconstruction
iZern = 1;
xi = 1: Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
wMin = -0.25;
wMax = 0.25;
waveFrontForShow = nan(Sx, Sx);
c = coeffs(iZern,:);
c1 = coeffsEst1(iZern,:);
c2 = coeffsEst2(iZern,:);
F1 = figure; % show wavefronts
figure(F1), subplot(1,3,1);
waveFrontForShow(idx) = create_wavefront(pIn,c,r(idx),theta(idx));
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
caxis([wMin wMax])
title('Input Wavefront');
figure(F1), subplot(1,3,2);
waveFrontForShow(idx) = create_wavefront(pIn,c1,r(idx),theta(idx));
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
caxis([wMin wMax])
title('Direct Recon');
figure(F1), subplot(1,3,3);
waveFrontForShow(idx) = create_wavefront(pIn,c2,r(idx),theta(idx));
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
caxis([wMin wMax])
title('Bootstrapping Recon');
savefig([fileFolderOut, fileNameOut, '_WF.fig']);



