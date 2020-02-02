% Test phase retrieval on experimental phase diversity images

clear all;
close all;
tic;
fileFolderIn = '.\DataForTest\';
fileFolderOut = '.\DataForTest\Output\';
fileImgIn = [fileFolderIn,'img.tif'];
fileZernIn = [fileFolderIn,'zern.dat'];
imgNum = 3; % 2 or 3
repNum = 1;
zernCoeffOrder = 21; % 4th (15), 5th(21), 6th(28) and 7th(36) orders
flagExcludeTilt = 1; % exclude tilts
flagGPU = 1;
if(flagExcludeTilt==1)
    pIn = 4:zernCoeffOrder; % 1: piston; 2:tilt X; 3: tilt Y;
else
    pIn = 2:zernCoeffOrder;
end
zernNum = length(pIn);
iteration = 5; % note: more zernike orders --> more iterations? *******
gamma = 1e-14;
zernSteps = [15 pIn(zernNum)] - pIn(1) + 1; % upboosting steps,e.g., 4th (15), 5th(21), 6th(28) and 7th(36) orders

pixelSize = 0.096; % um
lambda = 0.550; % um
NA = 1.2;
% create output folder
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
filePreImgs = [fileFolderOut 'imgs_pre.tif']; % image: pre rotation and cropping
fileImgs = [fileFolderOut 'imgs.tif']; % image

disp('...Preprocessing images...');
% % % input images
rotAng = 72;
cropSize = 512;
bgValue = 1000;
imgsRaw = double(ReadTifStack(fileImgSample));
[Sx1, Sy1, rawNum] = size(imgsRaw);
% average if acquistion repeated for each phase
if(repNum>1)
    aveNum = round(rawNum/repNum);
    imgsAve = zeros(Sx1,Sy1,rawNum);
    for i = 1:aveNum
        imgAve = zeros(Sx1,Sy1);
        for j = 1:repNum
            iSlice = (i-1)*repNum + j;
            if(iSlice <=rawNum)
                imgAve = imgAve + imgsRaw(:,:,(i-1)*repNum + j);
            end
        end
        imgsAve(:,:,i) = imgAve;
    end
    imgsAve = imgsAve/repNum;
else
    aveNum = rawNum;
    imgsAve = imgsRaw;
end
WriteTifStack(imgsAve, filePreImgs, 32);
% rotation and crop
imgs = zeros(cropSize,cropSize,imgNum);
for i = 1:imgNum
    img = imrotate(imgs(:,:,i),rotAng,'bilinear');
    [Sox,Soy] = round((Size(img)-cropSize)/2+1);
    imgs(:,:,i) = img(Sox:Sox+cropSize-1,Soy:Soy+cropSize-1);
end
imgs = imgs - bgValue;
imgs = max(imgs,0.01);
WriteTifStack(imgs, [fileFolderOut 'Image_phasediversity.tif'], 32);

% % % phases and zernike coefficients
coeffsRaw = importdata('test_zern.txt');
[phaseNum, zernRawNum] = size(coeffsRaw);
if(phaseNum<imgNum)
    error('phase number is less than image number');
end
if((zernRawNum + 1) < zernCoeffOrder)
    error('raw zernike number is less than configured zernike number');
end
% Imagine Opitc convetion to Fringe convention
% 1st: tilt <--> 1st: piston
if(flagExcludeTilt==1)
    pStart = 3; % 1:tilt X; 2: tilt Y; 3: defocus;
else
    pStart = 1;
end
% unknown aberration:
coeffsInitial = coeffsRaw(1,pStart:zernCoeffOrder);
% phase diversity:
coeffs_delta = coeffsRaw(2:imgNum,pStart:zernCoeffOrder);
coeffs_all = [coeffsInitial';coeffs_delta';]';

cTime1 = toc;
disp(['... ... time cost: ', num2str(cTime1)]);
% reconstruct zernike coefficients: estimate the unknown aberration
disp('...Reconstructing Zernike coefficients...');
devNum = 1;
gpuDevice(devNum);
[cEstimate, imgEstimate, rePar] = recon_zern(imgs, pIn, coeffs_delta, gamma, iteration, zernSteps, pixelSize, lambda, NA, flagGPU);
cTime2 = toc;
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
save([fileFolderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime2)]);

% calculate PSFs and waveFront based on input zernike coefficients
img0 = squeeze(imgs(:,:,1));
[~, PSFs, waveFronts] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, lambda, NA, 'none');
WriteTifStack(PSFs, [fileFolderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(imgEstimate, [fileFolderOut 'Image_estimated.tif'], 32);
[Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront = zeros(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront, [fileFolderOut 'Wavefront_estimated.tif'], 32);
% export zernike coefficient to txt file: imagine optic convention
% 1st: piston <--> 1st: tilt
coeffsOut = zeros(1,zernRawNum);
coeffsOut(1,pStart:zernCoeffOrder) = cEstimate;
fileID = fopen('zern_estimated.txt','w');
fprintf(fileID,'%f\t',coeffsOut);
fclose(fileID);
fileID = fopen('zern_estimated_negative.txt','w');
fprintf(fileID,'%f\t',-coeffsOut);
fclose(fileID);

% check results
xi = 1: Sx;
a = 20; % size of PSF images for show
F1 = figure; % input images and phases
figure(F1), subplot(imgNum,3,1);
pcolor(xi,xi,waveFronts(:,:,1)), shading interp
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
pcolor(xi,xi,waveFronts(:,:,2)), shading interp
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
    pcolor(xi,xi,waveFronts(:,:,3)), shading interp
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
pcolor(xi,xi,waveFronts(:,:,1)), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront: aberrated');

figure(F2); subplot(2,2,2);
img = img0/max(img0(:));
imshow(img,[]),colorbar;
title('Image: aberrated');

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
xlabel('Zernike Coeff Order');
ylabel('Zernike Coeff Magnitude');
set(gca,'FontSize', 14);
title('Zernike coefficients');
savefig([fileFolderOut 'coeff.fig']);

