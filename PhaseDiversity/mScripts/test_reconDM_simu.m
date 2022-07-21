% A simple example: test phase diveristy for DM calibration based on synthetic images.
%   1) generate diveristy images series;
%   2) reconstruct the initial and delta slope wavefront (Zernike coefficients);
%   3) show results.
%
% Current setup has been verified for the test dataset 'imSample1': 
% (import 'InputParameters.mat', and commet Line 32) 
% or users may need to cutomize their input aberrations and parameters

% By: Min Guo
% Last update: Mar. 24, 2022

% % *** default settings:
clear all;
% close all;
tStart = tic;
flagGPU = 1;
pixelSize = 0.08; % um
lambda = 0.532; % um
NA = 1.2;
pIn = 3:20; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);
nType = 'none';
SNR = 20;

% % ****** customize wavefronts here ****************
% % Unknown initail wavefronts: to be estimated 
% random zernike coefficients
aValue = 0.2; % um: small value
zernType = 'astig'; % 'defocus','astig','coma','trefoil','sphe','random1','random2'
coeffsInitial = gen_zern_coeffs(pIn,aValue,zernType);

% % Unknown delta slope wavefronts: to be estimated 
% random zernike coefficients
aValue = 1; % um
zernType = 'defocus'; % 'defocus','astig','coma','trefoil','sphe','random1','random2'
coeffsSlope = gen_zern_coeffs(pIn,aValue,zernType);

% diversity images based on actuator steps:
kStepNum = 3; % 
imgNum = 2*kStepNum + 1;
deltaVol = 0.5; % delta voltage on actuator

% % reconstruction settings
itLimit = 10; % max iteration number
gamma = 1e-10; % gamma range: 1e-14 ~ 1e-6
alpha = 0;
penalChioce = 1;

% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
% zernSteps = [14 20 27 pIn(zernNum)] - pIn(1) + 1; 
zernSteps = [14 pIn(zernNum)] - pIn(1) + 1; 

% % files and folders
fileFolderIn = '..\DataForTest\Simu\';
fileImgSample = [fileFolderIn 'imSample1.tif'];
img0 = single(ReadTifStack(fileImgSample));
[Sx, ~] = size(img0);
fileFolderOut = '..\DataForTest\SimuDM\TestResults\';
% fileFolderOut = '..\..\..\computAO\Data\DataForTest\Simu\TestResults\';
fileImgs = [fileFolderOut 'imgs.tif']; % image
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

% delta wavefronts:
% kSteps = [0 -kStepNum:-1 1:kStepNum]; % organized to match previous images format
kSteps = -kStepNum:kStepNum; % organized to match previous images format
coeffs_all = zeros(imgNum,zernNum, 'single'); % 1 initial wavefront + 2K delta
for i = 1:imgNum
    coeffs_all(i,:) = coeffsInitial + kSteps(i) * deltaVol * coeffsSlope;
end
for i = 2:imgNum
    coeffs_all(i,:) = coeffs_all(i,:) - coeffs_all(1,:);
end
% coeffs_all(1,:) = coeffs_all(1,:) + coeffsInitial; % replace with initial wavefront


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
% [cEstimate, bEstimate, imgEstimate, rePar] = recon_zern_DM(imgs, pIn, deltaVol, ...
%     gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);
% [cEstimate, bEstimate, imgEstimate, rePar] = recon_zern_DM(imgs, pIn, deltaVol, ...
%     gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce, 0, coeffsSlope);
[cEstimate, bEstimate, imgEstimate, rePar] = recon_zern_DM(imgs, pIn, deltaVol, ...
    gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce, 0, coeffsInitial);

cTime2 = toc(tStart);
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);

% calculate RMS and PV error
[~, ~, staParaInitial, ~] = coeffs2wavefront(pIn,coeffsSlope,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaEstimate, ~] = coeffs2wavefront(pIn,bEstimate,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaError, ~] = coeffs2wavefront(pIn,bEstimate-coeffsSlope,Sx,...
    pixelSize, lambda, NA, 0);
save([fileFolderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime2)]);
WriteTifStack(img0, [fileFolderOut 'Image_groundtruth.tif'], 32);
WriteTifStack(imgs, [fileFolderOut 'Image_phasediversity.tif'], 32);
WriteTifStack(PSFs, [fileFolderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(imgEstimate, [fileFolderOut 'Image_estimated.tif'], 32);
WriteTifStack(waveFronts, [fileFolderOut 'Wavefront_phasediversity.tif'], 32);
[Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront_c = zeros(Sx, Sy);
waveFront_c(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront_c, [fileFolderOut 'Wavefront_estimated_c.tif'], 32);
waveFront_b = zeros(Sx, Sy);
waveFront_b(idx) = create_wavefront(pIn,bEstimate,r(idx),theta(idx));
WriteTifStack(waveFront_b, [fileFolderOut 'Wavefront_estimated_b.tif'], 32);

% check results
for i = 2:imgNum
    waveFronts(:,:,i) = waveFronts(:,:,i) + waveFronts(:,:,1);
end
for i = [1:kStepNum kStepNum+2:imgNum]
    waveFronts(:,:,i) = waveFronts(:,:,i) - waveFronts(:,:,kStepNum+1);
end

xi = 1: Sx;
a = 20;
waveFrontForShow = nan(Sx, Sx);
F1 = figure; % input images and phases
figure(F1), subplot(3,3,1);
waveFrontTemp = waveFronts(:,:,kStepNum+1);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: initial');
figure(F1); subplot(3,3,2);
PSF = PSFs(:,:,kStepNum+1);
PSF = PSF/max(PSF(:));
PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF: initial');
figure(F1); subplot(3,3,3);
img = imgs(:,:,1);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image: initial');

figure(F1), subplot(3,3,4);
waveFrontTemp = waveFronts(:,:,kStepNum);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: step -1');
figure(F1); subplot(3,3,5);
PSF = PSFs(:,:,kStepNum);
PSF = PSF/max(PSF(:));
PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF: step -1');
figure(F1); subplot(3,3,6);
img = imgs(:,:,2);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image: step -1');

figure(F1), subplot(3,3,7);
waveFrontTemp = waveFronts(:,:,kStepNum+2);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: step 1');
figure(F1); subplot(3,3,8);
PSF = PSFs(:,:,kStepNum+2);
PSF = PSF/max(PSF(:));
PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF: step 1');
figure(F1); subplot(3,3,9);
img = imgs(:,:,3);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image: step 1');
savefig([fileFolderOut 'input.fig']);

waveFront_original = waveFronts(:,:,kStepNum+1);
wMin_c = min(waveFront_original(:));
wMax_c = max(waveFront_original(:));
F2 = figure; subplot(2,3,1);
waveFrontTemp = waveFront_original;
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin_c wMax_c])
xlabel(Fc,'\mum');
title('Wavefront: ground truth c');

figure(F2); subplot(2,3,4);
waveFrontTemp = waveFront_c;
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin_c wMax_c])
xlabel(Fc,'\mum');
title('Wavefront: estimated c');

waveFront = zeros(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,coeffsSlope,r(idx),theta(idx));
wMin_b = min(waveFront(:));
wMax_b = max(waveFront(:));
figure(F2); subplot(2,3,2);
waveFrontTemp = waveFront;
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin_b wMax_b])
xlabel(Fc,'\mum');
title('Wavefront: ground truth b');

figure(F2); subplot(2,3,5);
waveFrontTemp = waveFront_b;
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
caxis([wMin_b wMax_b])
xlabel(Fc,'\mum');
title('Wavefront: estimated b');

figure(F2); subplot(2,3,3);
img = img0/max(img0(:));
imshow(img,[]),colorbar;
title('Image:ground truth');

figure(F2); subplot(2,3,6);
img = imgEstimate/max(imgEstimate(:));
imshow(img,[]),colorbar;
title('Image: estimated');
savefig([fileFolderOut 'retrieval.fig']);

F3 = figure; % plot Zernike coefficients
subplot(1,2,1);
plot(pIn,coeffsInitial,pIn,cEstimate,'LineWidth',2); 
legend( 'groundtruth','estimated');
xlabel('Zernike Coeff Index');
ylabel('Zernike Coeff Magnitude');
set(gca,'FontSize', 14);
title('Zernike coefficients c');

figure(F3); subplot(1,2,2); % plot Zernike coefficients
plot(pIn,coeffsSlope,pIn,bEstimate,'LineWidth',2); 
legend( 'groundtruth','estimated');
xlabel('Zernike Coeff Index');
ylabel('Zernike Coeff Magnitude');
set(gca,'FontSize', 14);
title('Zernike coefficients b');
savefig([fileFolderOut 'coeffs.fig']);

