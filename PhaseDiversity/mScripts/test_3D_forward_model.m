% test 3D forward model and phase retrieval on simulated images

% By: Min Guo
% Apr. 10, 2020;

% clear all;
% close all;
tic;

flagGPU = 0;

pIn = 4:21; % 1: piston; 2:tilt X; 3: tilt Y;
zernNum = length(pIn);
iteration = 5; % note: more zernike orders --> more iterations? *******
gamma = 1e-14;
% zernSteps = pIn(zernNum) - pIn(1) + 1;
% zernSteps = [15 21 28 pIn(zernNum)] - pIn(1) + 1; % upboosting steps,e.g., 4th (15), 5th(21), 6th(28) and 7th(36) orders
zernSteps = [15 pIn(zernNum)] - pIn(1) + 1; 
% random zernike coefficients
a = -0.05; b = 0.05;
% cRandom = a + (b-a).*rand(1,zernNum); 
% unknown aberration:
coeffsInitial = cRandom;
% phase diversity:
imgNum = 2; % 2 or 3
abType1 = 'ast'; % 'spherical' or 'defocus', or 'ast' (astigmatism)
abValue1 = 0.8; % um
% abType1 = 'defocus';
% abValue1 = 0.5; % um
abType2 = 'defocus';
abValue2 = -0.2; % um
% fileFolderIn = '..\DataForTest\Simu\';
% fileFolderOut = '..\DataForTest\Simu\TestResults\';
fileFolderIn = 'C:\Programs\computAO\Data_ImSample2\';
fileFolderOut = 'C:\Programs\computAO\Retrieval\forwardmodeltest\forward3D_ast0p8_repeat\';
% fileFolderOut = 'C:\Programs\computAO\Retrieval\forwardmodeltest\forward3D_de0p5_repeat\';
nType = 'none';
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
fileImgSample = [fileFolderIn 'imSample2_3D.tif'];
fileImgs = [fileFolderOut 'imgs.tif']; % image

pixelSize = 0.1625; % um
lambda = 0.510; % um
NA = 0.8;
% Known aberration: image 1
if (strcmp(abType1, 'defocus'))
    abIdx = 4; % 4th: defocus, 9th: spherical
elseif (strcmp(abType1, 'ast'))
    abIdx = 5; % 5th: astigmatism
elseif (strcmp(abType1, 'spherical'))
    abIdx = 9; % 4th: defocus, 9th: spherical
else
    error('No matched aberration type, please use defocus or spherical');
end
coeffs_delta1 = zeros(1,zernNum);
coeffs_delta1(abIdx- pIn(1) + 1) = abValue1; % um

coeffs_delta = coeffs_delta1;
% Known aberration: image 2
if (imgNum==3)
    if (strcmp(abType2, 'defocus'))
        abIdx = 4; % 4th: defocus, 9th: spherical
    elseif (strcmp(abType2, 'ast'))
        abIdx = 5; % 5th: astigmatism
    elseif (strcmp(abType2, 'spherical'))
        abIdx = 9; % 4th: defocus, 9th: spherical
    else
        error('No matched aberration type, please use defocus or spherical');
    end
    coeffs_delta2 = zeros(1, zernNum);
    coeffs_delta2(abIdx - pIn(1) + 1) = abValue2; % um
    coeffs_delta = [coeffs_delta;coeffs_delta2;];
end
coeffs_all = [coeffsInitial;coeffs_delta;];

% generate phase diversity images;
disp('...Generating simulated images...');
img0 = double(ReadTifStack(fileImgSample));
[Sx, Sy, Sz] = size(img0);
zSliceIdx = 80;
RI = 1.33;
zRatio = 5;
if(imgNum ==2)
    imgs = zeros(Sx,Sy,2);
else
    imgs = zeros(Sx,Sy,3);
end
[img3D1, PSF3D1] = gen_simu_3Dimage(img0, pIn, coeffsInitial, pixelSize, ...
    lambda, NA, pixelSize*zRatio, RI, nType);
WriteTifStack(img3D1,[fileFolderOut 'img3D1.tif'],32);
WriteTifStack(PSF3D1,[fileFolderOut 'PSF3D1.tif'],32);
imgs(:,:,1) = img3D1(:,:,zSliceIdx);

[img3D2, PSF3D2] = gen_simu_3Dimage(img0, pIn, coeffsInitial+coeffs_delta1, pixelSize, ...
    lambda, NA, pixelSize*zRatio, RI, nType);
WriteTifStack(img3D2,[fileFolderOut 'img3D2.tif'],32);
WriteTifStack(PSF3D2,[fileFolderOut 'PSF3D2.tif'],32);
imgs(:,:,2) = img3D2(:,:,zSliceIdx);
if(imgNum==3)
    [img3D3, PSF3D3] = gen_simu_3Dimage(img0, pIn, coeffsInitial+coeffs_delta1, pixelSize, ...
        lambda, NA, pixelSize*zRatio, RI, nType);
    WriteTifStack(img3D3,[fileFolderOut 'img3D3.tif'],32);
    WriteTifStack(PSF3D3,[fileFolderOut 'PSF3D3.tif'],32);
    imgs(:,:,3) = img3D3(:,:,zSliceIdx);
end
cTime1 = toc;
disp(['... ... time cost: ', num2str(cTime1)]);
[imgs2, PSFs, waveFronts] = gen_simu_images(img0(:,:,zSliceIdx), pIn, coeffs_all, pixelSize, lambda, NA, nType);

% reconstruct zernike coefficients: estimate the unknown aberration
disp('...Reconstructing Zernike coefficients...');
if(flagGPU)
    devNum = 1;
    gpuDevice(devNum);
end
[cEstimate, imgEstimate, rePar] = recon_zern(imgs, pIn, coeffs_delta, gamma, iteration, zernSteps, pixelSize, lambda, NA, flagGPU);
cTime2 = toc;
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
save([fileFolderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime2)]);
WriteTifStack(img0, [fileFolderOut 'Image_groundtruth.tif'], 32);
WriteTifStack(imgs, [fileFolderOut 'Image_phasediversity.tif'], 32);
WriteTifStack(PSFs, [fileFolderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(imgEstimate, [fileFolderOut 'Image_estimated.tif'], 32);
WriteTifStack(waveFronts, [fileFolderOut 'Wavefront_phasediversity.tif'], 32);
% [Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront = zeros(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront, [fileFolderOut 'Wavefront_estimated.tif'], 32);

% check results
xi = 1: Sx;
a = 20;
F1 = figure; % input images and phases
figure(F1), subplot(imgNum,3,1);
pcolor(xi,xi,waveFronts(:,:,1)), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: aberrated');
figure(F1); subplot(imgNum,3,2);
PSF0 = PSFs(:,:,1);
PSF0 = PSF0/max(PSF0(:));
PSF = PSF0(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF:aberrated');
figure(F1); subplot(imgNum,3,3);
img = imgs(:,:,1);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image:aberrated');

figure(F1), subplot(imgNum,3,4);
pcolor(xi,xi,waveFronts(:,:,2)), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront:phase1');
figure(F1); subplot(imgNum,3,5);
PSF0 = PSFs(:,:,2);
PSF0 = PSF0/max(PSF0(:));
PSF = PSF0(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
imshow(PSF,[]),colorbar;
title('PSF:phase1');
figure(F1); subplot(imgNum,3,6);
img = imgs(:,:,2);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image:phase1');

if(imgNum ==3)
    figure(F1), subplot(imgNum,3,7);
    pcolor(xi,xi,waveFronts(:,:,3)), shading interp
    axis square, Fc = colorbar;
    xlabel(Fc,'\mum');
    title('Wavefront:phase2');
    figure(F1); subplot(imgNum,3,8);
    PSF0 = PSFs(:,:,3);
    PSF0 = PSF0/max(PSF0(:));
    PSF = PSF0(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
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
pcolor(xi,xi,waveFronts(:,:,1)), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: aberrated');

figure(F2); subplot(2,2,2);
img = img0(:,:,zSliceIdx);
img = img/max(img(:));
imshow(img,[]),colorbar;
title('Image:ground truth');

figure(F2); subplot(2,2,3);
pcolor(xi,xi,waveFront), shading interp
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

