% Explores the effect from diversity phase noise

% Jul. 29, 2020

% clear all;
% close all;
tic;
flagGPU = 0;

pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);
itLimit = 10; % bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
gamma = 1e-14;
% alpha = 1e4;
alpha = 0;
nType = 'none';
% nType = 'both';
SNR = 20;
deltaNoiseIdx = 1;
deltaNoiseValue = 0;
folderName = 'NoiseFree_test_imSample1';

% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
zernSteps = [14 20 27 pIn(zernNum)] - pIn(1) + 1; 

% random zernike coefficients
aValue = 0.1;
% % % uniform weight for all indices
% zernType = 'random2';
% % % higher weight for lower order indices
zernType = 'random2';
% coeffsInitial = gen_zern_coeffs(pIn,aValue,zernType);

% phase diversity:
imgNum = 2; % 2 or 3
abType1 = 'defocus'; % 'astig' or 'defocus'
abType2 = 'astig';
abValue1 = 0.6; % um
abValue2 = -0.2; % um
fileFolderIn = '..\DataForTest\Simu\';
fileFolderOut = ['..\..\..\computAO\Data\exploreDeltaNoise\', folderName,'\'];
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
fileImgSample = [fileFolderIn 'imSample1.tif'];
fileImgs = [fileFolderOut 'imgs.tif']; % image

pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
% Known aberration: image 1
coeffs_delta1 = gen_zern_coeffs(pIn,abValue1,abType1);

coeffs_delta = coeffs_delta1;
% Known aberration: image 2
if (imgNum==3)
    coeffs_delta2 = gen_zern_coeffs(pIn,abValue2,abType2);
    coeffs_delta = [coeffs_delta;coeffs_delta2;];
end
coeffs_all = [coeffsInitial;coeffs_delta;];

pNoise = 1:max(pIn);
cLengthNoise = length(pNoise);
coeffs_deltaNoise = zeros(size(coeffs_delta,1),cLengthNoise);
coeffs_deltaNoise(:,pIn(1):cLengthNoise) = coeffs_delta;
coeffs_deltaNoise(1,deltaNoiseIdx) = coeffs_deltaNoise(1,deltaNoiseIdx) + deltaNoiseValue; % um
coeffsInitialNoise(1,pIn(1):cLengthNoise) = coeffsInitial;
coeffs_allNoise = [coeffsInitialNoise;coeffs_deltaNoise;];

% generate phase diversity images;
disp('...Generating simulated images...');
img0 = double(ReadTifStack(fileImgSample));
% [imgs, PSFs, waveFronts] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, lambda, NA, nType);
[imgs, PSFs, waveFronts] = gen_simu_images(img0, pNoise, coeffs_allNoise, pixelSize, lambda, NA, nType,SNR);
cTime1 = toc;
disp(['... ... time cost: ', num2str(cTime1)]);
% reconstruct zernike coefficients: estimate the unknown aberration
disp('...Reconstructing Zernike coefficients...');
if(flagGPU)
    devNum = 1;
    gpuDevice(devNum);
end
[cEstimate, imgEstimate, rePar] = recon_zern(imgs, pIn, coeffs_delta, ...
    gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU);
cTime2 = toc;
% RMS and PV
[~, ~, staParaInitial, ~] = coeffs2wavefront(pIn,coeffsInitial,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaEstimate, ~] = coeffs2wavefront(pIn,cEstimate,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaError, ~] = coeffs2wavefront(pIn,cEstimate-coeffsInitial,Sx,...
    pixelSize, lambda, NA, 0);
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
save([fileFolderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime2)]);
WriteTifStack(img0, [fileFolderOut 'Image_groundtruth.tif'], 32);
WriteTifStack(imgs, [fileFolderOut 'Image_phasediversity.tif'], 32);
WriteTifStack(PSFs, [fileFolderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(waveFronts, [fileFolderOut 'Wavefront_phasediversity.tif'], 32);
WriteTifStack(imgEstimate, [fileFolderOut 'Image_estimated.tif'], 32);
[Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront = zeros(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront, [fileFolderOut 'Wavefront_estimated.tif'], 32);

% % % check results
% waveFront_original = waveFronts(:,:,1);
% wMin = min(waveFront_original(:));
% wMax = max(waveFront_original(:));
wMin = -0.25;
wMax = 0.25;
waveFrontForShow = nan(Sx, Sy);
xi = 1: Sx;
a = 20;
F1 = figure; % input images and phases
figure(F1), subplot(imgNum,3,1);
waveFrontTemp = waveFronts(:,:,1);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
caxis([wMin wMax])
title('Wavefront: aberrated');
figure(F1); subplot(imgNum,3,2);
PSF = PSFs(:,:,1);
PSF = PSF/max(PSF(:));
PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF:aberrated');
figure(F1); subplot(imgNum,3,3);
img = imgs(:,:,1);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image:aberrated');

figure(F1), subplot(imgNum,3,4);
waveFrontTemp = waveFronts(:,:,2);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront:phase1');
figure(F1); subplot(imgNum,3,5);
PSF = PSFs(:,:,2);
PSF = PSF/max(PSF(:));
PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF:phase1');
figure(F1); subplot(imgNum,3,6);
img = imgs(:,:,2);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image:phase1');

if(imgNum ==3)
    figure(F1), subplot(imgNum,3,7);
    waveFrontTemp = waveFronts(:,:,3);
    waveFrontForShow(idx) = waveFrontTemp(idx);
    pcolor(xi,xi,waveFrontForShow), shading interp
    axis square, Fc = colorbar;
    xlabel(Fc,'\mum');
    title('Wavefront:phase2');
    figure(F1); subplot(imgNum,3,8);
    PSF = PSFs(:,:,3);
    PSF = PSF/max(PSF(:));
    PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
    imshow(PSF,[]),colorbar;
    title('PSF:phase2');
    figure(F1); subplot(imgNum,3,9);
    img = imgs(:,:,3);
    img = img/max(img(:));
    imshow(img,[]),colorbar;
    title('Image:phase2');
end
savefig([fileFolderOut 'input.fig']);



F2 = figure; subplot(2,2,1);
waveFrontTemp = waveFronts(:,:,1);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: aberrated');

figure(F2); subplot(2,2,2);
img = img0/max(img0(:));
imshow(img,[]),colorbar;
title('Image:ground truth');

figure(F2); subplot(2,2,3);
waveFrontForShow(idx) = waveFront(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: estimated');

figure(F2); subplot(2,2,4);
img = imgEstimate/max(imgEstimate(:));
imshow(img,[]),colorbar;
title('Image:estimated');
savefig([fileFolderOut 'retrieval.fig']);

F3 = figure; % plot Zernike coefficients
plot(pIn,coeffsInitial,pIn,cEstimate,'LineWidth',2); 
legend( 'groundtruth','estimated');
xlabel('Zernike Index');
ylabel('Coeff Value');
set(gca,'FontSize', 14);
title('Zernike coefficients');
savefig([fileFolderOut 'coeff.fig']);
