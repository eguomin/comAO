% A simple example: test phase diveristy based on synthetic images.
%   1) generate aberrated image (and diveristy images);
%   2) reconstruct the aberrated wavefront (Zernike coefficients);
%   3) show results.
%
% Current setup has been verified for the test dataset 'imSample1': 
% (import 'InputParameters.mat', and commet Line 32) 
% or users may need to cutomize their input aberrations and parameters

% By: Min Guo
% Last update: Mar. 11, 2022

% % *** default settings:
% clear all;
% close all;
tStart = tic;
flagGPU = 0;
pixelSize = 0.08; % um
lambda = 0.532; % um
NA = 1.2;
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);
nType = 'none';
SNR = 20;

% % ****** customize aberrations here ****************
% % Unknown aberrations: to be estimated 
% random zernike coefficients
aValue = 0.2;
zernType = 'astig'; % 'defocus','astig','coma','trefoil','sphe','random1','random2'
coeffsInitial = gen_zern_coeffs(pIn,aValue,zernType);


% diversity images:
imgNum = 3; % 2 or 3
abType1 = 'defocus'; % 'astig' or 'defocus'
abValue1 = 0.5; % um
abType2 = 'defocus'; % astig at 45 deg
abValue2 = -0.5; % um
zernIdx = 4; 
zernValue = abValue1;

% % reconstruction settings
itLimit = 10; % max iteration number
gamma = 1e-10; % gamma range: 1e-14 ~ 1e-6
alpha = 0;
penalChioce = 1;

% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
zernSteps = [14 20 27 pIn(zernNum)] - pIn(1) + 1; 

% % files and folders
fileFolderIn = '..\DataForTest\Simu\';
fileImgSample = [fileFolderIn 'imSample1.tif'];
img0 = single(ReadTifStack(fileImgSample));
[Sx, ~] = size(img0);
fileFolderOut = '..\DataForTest\Simu\TestResults2\';
% fileFolderOut = '..\..\..\computAO\Data\DataForTest\Simu\TestResults\';
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
[cEstimate, imgEstimate, rePar] = recon_zern(imgs, pIn, coeffs_delta, ...
    gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);
cTime2 = toc(tStart);
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);

% calculate RMS and PV error
[~, ~, staParaInitial, ~] = coeffs2wavefront(pIn,coeffsInitial,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaEstimate, ~] = coeffs2wavefront(pIn,cEstimate,Sx,...
    pixelSize, lambda, NA, 0);
[~, ~, staParaError, ~] = coeffs2wavefront(pIn,cEstimate-coeffsInitial,Sx,...
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
waveFront = zeros(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront, [fileFolderOut 'Wavefront_estimated.tif'], 32);

% check results
xi = 1: Sx;
a = 20;
waveFrontForShow = nan(Sx, Sx);
F1 = figure; % input images and phases
figure(F1), subplot(imgNum,3,1);
waveFrontTemp = waveFronts(:,:,1);
waveFrontForShow(idx) = waveFrontTemp(idx);
pcolor(xi,xi,waveFrontForShow), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
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

waveFront_original = waveFronts(:,:,1);
wMin = min(waveFront_original(:));
wMax = max(waveFront_original(:));
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
waveFrontTemp = waveFront;
waveFrontForShow(idx) = waveFrontTemp(idx);
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
xlabel('Zernike Coeff Index');
ylabel('Zernike Coeff Magnitude');
set(gca,'FontSize', 14);
title('Zernike coefficients');
savefig([fileFolderOut 'coeff.fig']);

