function [coeffsOut, cTime] = recon_zern_fileIO(fileFolderOut,fileFolderIn, fileName, imgNum, repNum,...
    zernCoeffOrder, iteration, gamma, cropSize, bgValue,flagShowInput, flagShowRecon, flagGPU)
% phase retrieval based on input files  of phase diversity images
%  % Output
% coeffsOut: reconstructed Zernike coefficients, to be applied to DM to
% cancel aberration;
% cTime: computational time;
% % Input
% fileFolderOut: strings of the output path;
% fileFolderIn: strings of the input path;
% fileName: strings of the file name;
% imgNum: number of phase diversity images (including raw image);
% repNum: number of images for each phase diversity;
% zernCoeffOrder: Zernike coefficient index to be retrieved;
% iteration: iteration number for the reconstruction;
% gamma: tuning parameter, 1e-5~1e-10;
% cropSize: image size for calculation;
% bgValue: background of the input image;
% flagShowInput: show image/phase of input;
% flagShowRecon: show reconstruction result;
% flagGPU: use GPU or not;
% by Min Guo
% Feb. 9, 2020;

% clear all;
% close all;
warning('off','all')
tic;
% fileFolderIn = 'D:\Data\20200207_PR\';
% fileName = 'beads_2';
% fileFolderOut = [fileFolderIn,fileName, '\'];
fileImgIn = [fileFolderIn, fileName,'.tif'];
fileZernIn = [fileFolderIn,fileName,'_coeff.txt'];
% imgNum = 5; % 2 or 3
% fileFolderOut = [fileFolderIn,fileName, num2str(imgNum),'\'];
% repNum = 2;
% zernCoeffOrder = 15; % 4th (15), 5th(21), 6th(28) and 7th(36) orders
flagExcludeTilt = 1; % exclude tilts
flagExcludeDefocus = 0; % exclude defocus
rotAng = 70;
% flagShowInput = 0;
% flagShowRecon = 1;
flagSmoothEdge = 1;

% flagGPU = 1;
if(flagExcludeTilt==1)
    pIn = 4:zernCoeffOrder; % 1: piston; 2:tilt X; 3: tilt Y;
else
    pIn = 2:zernCoeffOrder;
end
if(flagExcludeDefocus==1)
    pIn = 5:zernCoeffOrder; % 1: piston; 2:tilt X; 3: tilt Y;
end

zernNum = length(pIn);
% iteration = it; % note: more zernike orders --> more iterations? *******
% gamma = 1e-6;

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

disp('...Preprocessing images...');
% % % input images
% rotAng = -73;
% cropSize = 384;
% bgValue = 290;
imgsRaw = double(ReadTifStack(fileImgIn));
[Sx1, Sy1, rawNum] = size(imgsRaw);
% average if acquistion repeated for each phase
if(repNum>1)
    aveNum = round(rawNum/repNum);
    imgsAve = zeros(Sx1,Sy1,aveNum);
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
if(aveNum<imgNum)
    error('image number in TIFF is less than image number');
end
imgs = zeros(cropSize,cropSize,imgNum);
for i = 1:imgNum
    imgIn = imgsAve(:,:,i);
    imgIn = flipud(imgIn);
    img = imrotate(imgIn,rotAng,'bilinear');
%     img = flipud(img);
    Soxy = round((size(img)-cropSize)/2+1);
    Sox = Soxy(1);
    Soy = Soxy(2);
    imgs(:,:,i) = img(Sox:Sox+cropSize-1,Soy:Soy+cropSize-1);
end
imgs = imgs - bgValue;
imgs = max(imgs,0.01);
% smooth the boundary of images
if(flagSmoothEdge==1)
    sKernel = fspecial('gaussian',[10 10],3); 
    for i = 1:imgNum
        img = imgs(:,:,i);
        img = edgetaper(img,sKernel);
        imgs(:,:,i) = img;
    end
end  
WriteTifStack(imgs, [fileFolderOut 'Image_phasediversity.tif'], 32);

% % % phases and zernike coefficients
coeffsRaw = -importdata(fileZernIn);
[phaseNum, zernRawNum] = size(coeffsRaw);
if(phaseNum<imgNum)
    error('phase number is less than image number');
end
if(zernRawNum < (zernCoeffOrder-1))
    error('raw zernike number is less than configured zernike number');
end
% coeffsRaw(2:3,:) = coeffsRaw(4:5,:);
coeffsSigns = ones(1,zernRawNum); % to correct the mismatch between MATLAB and HASO
coeffsSigns([2,5, 7, 10, 12, 15]) = -1; % need further update
% zernSigns = -zernSigns;
for i = 1:phaseNum
    coeffsRaw(i,:) = coeffsRaw(i,:).*coeffsSigns;
end
% Imagine Opitc convetion to Fringe convention
% 1st: tilt <--> 1st: piston
if(flagExcludeTilt==1)
    pStart = 3; % 1:tilt X; 2: tilt Y; 3: defocus;
else
    pStart = 1;
end
if(flagExcludeDefocus==1)
    pStart = 4;
end
pEnd = zernCoeffOrder-1;
% unknown aberration:
coeffsInitial = coeffsRaw(1,pStart:pEnd);
coeffsInitial = coeffsInitial';
% phase diversity:
coeffs_delta = coeffsRaw(2:imgNum,pStart:pEnd);
coeffs_delta = coeffs_delta';
coeffs_all = [coeffsInitial';coeffs_delta';]';

cTime1 = toc;
disp(['... ... time cost: ', num2str(cTime1)]);
% reconstruct zernike coefficients: estimate the unknown aberration
disp('...Reconstructing Zernike coefficients...');
devNum = 1;
% gpuDevice(devNum);
[cEstimate, imgEstimate, rePar] = recon_zern(imgs, pIn, coeffs_delta, gamma, iteration, zernSteps, pixelSize, lambda, NA, flagGPU);
cTime2 = toc;
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
save([fileFolderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime2)]);

% calculate PSFs and waveFront based on input zernike coefficients
img0 = squeeze(imgs(:,:,1));
[~, PSFs, waveFronts] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, lambda, NA, 'none');
WriteTifStack(PSFs, [fileFolderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(waveFronts, [fileFolderOut 'Wavefront_phasediversity.tif'], 32);
WriteTifStack(imgEstimate, [fileFolderOut 'Image_estimated.tif'], 32);
[Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront = zeros(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront, [fileFolderOut 'Wavefront_estimated.tif'], 32);
% export zernike coefficient to txt file: imagine optic convention
% 1st: piston <--> 1st: tilt
coeffsOut = zeros(1,zernRawNum);
coeffsOut(1,pStart:pEnd) = cEstimate;
coeffsOut = coeffsOut.*coeffsSigns;
fileID = fopen([fileFolderOut, 'zernCoeffs_estimated.txt'],'w');
fprintf(fileID,'%f\t',coeffsOut');
fclose(fileID);

% check results
xi = 1: Sx;

if(flagShowInput == 1)
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
title('Add:phase1');
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

if(imgNum >=3)
for i = 3:imgNum
    figure(F1), subplot(imgNum,3,3*(i-1)+1);
    pcolor(xi,xi,waveFronts(:,:,i)), shading interp
    axis square, Fc = colorbar;
    xlabel(Fc,'\mum');
    title(['Add:phase', num2str(i-1)]);
    figure(F1); subplot(imgNum,3,3*(i-1)+2);
    PSF = PSFs(:,:,i);
    PSF = PSF/max(PSF(:));
    PSF = PSF(round(Sx/2)-a:round(Sx/2)+a,round(Sy/2)-a:round(Sy/2)+a);
    imshow(PSF,[]),colorbar;
    title(['PSF:phase', num2str(i-1)]);
    figure(F1); subplot(imgNum,3,3*(i-1)+3);
    img = imgs(:,:,i);
    img = img/max(img(:));
    imshow(img,[]),colorbar;
    title(['Image:phase', num2str(i-1)]);
end
end
savefig([fileFolderOut 'input.fig']);
end
if(flagShowRecon == 1)
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
plot(pIn-1,coeffsInitial,pIn-1,cEstimate,'LineWidth',2); % imagine optic convention
legend( 'groundtruth','estimated');
xlabel('Zernike Coeff Order');
ylabel('Zernike Coeff Magnitude');
set(gca,'FontSize', 14);
title('Zernike coefficients');
savefig([fileFolderOut 'coeff.fig']);
end
cTime = toc;

