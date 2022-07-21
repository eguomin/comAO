% Test graduated optimization scheme based on simulated phase diversity images

% By: Min Guo
% Last update: Sep. 1, 2020

% clear all;
close all;
% % default settings:
tStart = tic;
flagGPU = 0;
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);
nType = 'both';
SNR = 20;

% % Unknown aberrations: to be estimated 
% random zernike coefficients
aValue = 0.15;
zernType = 'random2'; % 'defocus','astig','coma','trefoil','sphe','random1','random2'
coeffsInitial = gen_zern_coeffs(pIn,aValue,zernType);
% coeffsInitial = c;

% diversity images:
imgNum = 2; % 2 or 3
abType1 = 'astig'; % 'astig' or 'defocus'
abValue1 = 0.5; % um
abType2 = 'idx'; % astig at 45 deg
abValue2 = 0.5; % um
zernIdx = 3; 
zernValue = abValue1;

% % reconstruction settings
itLimit = 10; % max iteration number
gamma = 1e-10; % % gamma range: 1e-14 ~ 1e-6
alpha = 0;
penalChioce = 1;
zernSteps1 = pIn(zernNum) - pIn(1) + 1;
% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
zernSteps2 = [14 20 27 pIn(zernNum)] - pIn(1) + 1; 

% % files and folders
fileFolderIn = '..\DataForTest\Simu\';
fileImgSample = [fileFolderIn 'imSample1.tif'];
img0 = single(ReadTifStack(fileImgSample));
[Sx, ~] = size(img0);
fileFolderOut = '..\..\..\computAO\Data\bootstrapping_new\none_astig_0p5\';
fileImgs = [fileFolderOut 'imgs.tif']; % image
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

% Known aberration: for diverisity image 1
coeffs_delta1 = gen_zern_coeffs(pIn,abValue1,abType1);
coeffs_delta = coeffs_delta1;

% Known aberration: for diverisity image 2
if (imgNum==3)
    coeffs_delta2 = gen_zern_coeffs(pIn,abValue2,abType2,zernIdx,zernValue);
    coeffs_delta = [coeffs_delta;coeffs_delta2;];
end
coeffs_all = [coeffsInitial;coeffs_delta;];

% generate phase diversity images;
disp('...Generating simulated images...');
[imgs, PSFs, waveFronts] = gen_simu_images(img0, pIn, coeffs_all, ...
    pixelSize, lambda, NA, nType, SNR);
cTime1 = toc(tStart);
disp(['... ... time cost: ', num2str(cTime1)]);
% reconstruct zernike coefficients: estimate the unknown aberration
disp('...Reconstructing Zernike coefficients...');
if(flagGPU)
    devNum = 1;
    gpuDevice(devNum);
end
disp('...direct reconstruction...');
[cEstimate1, imgEstimate1, rePar1] = recon_zern(imgs, pIn, coeffs_delta, ...
    gamma, alpha, itLimit, zernSteps1, pixelSize, lambda, NA, flagGPU, penalChioce);
cTime2 = toc(tStart);
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);

disp('...bootstrapping reconstruction...');
[cEstimate2, imgEstimate2, rePar2] = recon_zern(imgs, pIn, coeffs_delta, ...
    gamma, alpha, itLimit, zernSteps2, pixelSize, lambda, NA, flagGPU, penalChioce);
cTime3 = toc(tStart);
disp(['... ... time cost: ', num2str(cTime3-cTime2)]);

[~, ~, staParaInitial, ~] = coeffs2wavefront(pIn,coeffsInitial,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaEstimate1, ~] = coeffs2wavefront(pIn,cEstimate1,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaError1, ~] = coeffs2wavefront(pIn,cEstimate1-coeffsInitial,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaEstimate2, ~] = coeffs2wavefront(pIn,cEstimate2,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaError2, ~] = coeffs2wavefront(pIn,cEstimate2-coeffsInitial,Sx,...
    pixelSize, lambda, NA, 0);
save([fileFolderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime3)]);
WriteTifStack(img0, [fileFolderOut 'Image_groundtruth.tif'], 32);
WriteTifStack(imgs, [fileFolderOut 'Image_phasediversity.tif'], 32);
WriteTifStack(PSFs, [fileFolderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(waveFronts, [fileFolderOut 'Wavefront_phasediversity.tif'], 32);
WriteTifStack(imgEstimate1, [fileFolderOut 'Image_estimated_direct.tif'], 32);
WriteTifStack(imgEstimate2, [fileFolderOut 'Image_estimated_upboosting.tif'], 32);
[Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront1 = zeros(Sx, Sy);
waveFront1(idx) = create_wavefront(pIn,cEstimate1,r(idx),theta(idx));
WriteTifStack(waveFront1, [fileFolderOut 'Wavefront_estimated_direct.tif'], 32);
waveFront2 = zeros(Sx, Sy);
waveFront2(idx) = create_wavefront(pIn,cEstimate2,r(idx),theta(idx));
WriteTifStack(waveFront2, [fileFolderOut 'Wavefront_estimated_upboosting.tif'], 32);

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

F2 = figure; subplot(2,3,1);
waveFrontTemp = waveFronts(:,:,1);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: aberrated');

figure(F2); subplot(2,3,4);
img = img0/max(img0(:));
imshow(img,[]),colorbar;
title('Image:ground truth');

figure(F2); subplot(2,3,2);
waveFrontForShow(idx) = waveFront1(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: direct');

figure(F2); subplot(2,3,5);
img = imgEstimate1/max(imgEstimate1(:));
imshow(img,[]),colorbar;
title('Image:direct');

figure(F2); subplot(2,3,3);
waveFrontForShow(idx) = waveFront2(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: bootstrapping');

figure(F2); subplot(2,3,6);
img = imgEstimate2/max(imgEstimate2(:));
imshow(img,[]),colorbar;
title('Image:upboosting');
savefig([fileFolderOut 'retrieval.fig']);

F3 = figure; % plot Zernike coefficients
plot(pIn,coeffsInitial,'LineWidth',4); 
hold on;
plot(pIn,cEstimate1,pIn,cEstimate2,'LineWidth',2); 
legend( 'groundtruth','direct','bootstrapping');
xlabel('Zernike Index');
ylabel('Coeff Value');
set(gca,'FontSize', 14);
% title('Zernike coefficients');
savefig([fileFolderOut 'coeff.fig']);

F4 = figure; subplot(2,3,1);
waveFrontTemp = waveFronts(:,:,1);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: aberrated');
figure(F4); subplot(2,3,2);
waveFrontForShow(idx) = waveFront1(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: direct');
figure(F4); subplot(2,3,3);
waveFrontForShow(idx) = waveFront2(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: bootstrapping');
figure(F4); subplot(2,3,4);
cTemp = rePar2.cEstimate(:,1);
waveFrontForShow(idx) = create_wavefront(pIn,cTemp(:),r(idx),theta(idx));
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: step1');
figure(F4); subplot(2,3,5);
cTemp = rePar2.cEstimate(:,2);
waveFrontForShow(idx) = create_wavefront(pIn,cTemp(:),r(idx),theta(idx));
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: step2');
figure(F4); subplot(2,3,6);
cTemp = rePar2.cEstimate(:,3);
waveFrontForShow(idx) = create_wavefront(pIn,cTemp(:),r(idx),theta(idx));
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: step3');
savefig([fileFolderOut 'wavefronts_bootstrapping.fig']);

F5 = figure;
pIn1 = pIn(1:zernSteps2(1));
cTemp1 = rePar2.cEstimate(1:zernSteps2(1),1);
pIn2 = pIn(1:zernSteps2(2));
cTemp2 = rePar2.cEstimate(1:zernSteps2(2),2);
pIn3 = pIn(1:zernSteps2(3));
cTemp3 = rePar2.cEstimate(1:zernSteps2(3),3);
plot(pIn,coeffsInitial,pIn1,cTemp1(:),pIn2,cTemp2(:),...
    pIn3,cTemp3(:),pIn,cEstimate2,'LineWidth',2); 
legend( 'groundtruth','bootstrapping-step1','bootstrapping-step2',...
    'bootstrapping-step3','bootstrapping-final');
xlabel('Zernike Index');
ylabel('Coeff Value');
set(gca,'FontSize', 14);
title('Zernike coefficients');
savefig([fileFolderOut 'coeff_bootstrapping.fig']);