function [coeffsOut, cRMSE, cTime] = recon_zern_fileIO(folderOut,folderIn, fileName, ...
    imgNum, repNum, zernIdxLimit, itLimit, gamma, alpha, cropSize, bgValue,flagShowInput,...
    flagShowRecon, flagGPU,idxType,penalChioce)
% phase retrieval based on input files  of phase diversity images
%  % Output
% coeffsOut: reconstructed Zernike coefficients, to be applied to DM to
% cancel aberration;
% cRMSE: root mean squared error;
% cTime: computational time;
% % Input
% folderOut: strings of the output path;
% folderIn: strings of the input path;
% fileName: strings of the file name;
% imgNum: number of phase diversity images (including raw image);
% repNum: number of images for each phase diversity;
% zernIdxLimit: Zernike coefficient index to be retrieved;
% itLimit: maximum iteration number for the reconstruction;
% gamma: tuning parameter for image penalty term, 1e-5~1e-10;
% alpha: tuning parameter for phase penalty term, 0~1e8;
% cropSize: image size for calculation;
% bgValue: background of the input image;
% flagShowInput: show image/phase of input;
% flagShowRecon: show reconstruction result;
% flagGPU: use GPU or not;

% by Min Guo
% Feb. 9, 2020;
% Modification (Feb. 13, 2020):
%   add the calculation of RMSE value
% Modification: July 1, 2020
%   penalChioce: penalty term to image, 
%       1 for L2 norm of image and 2 for L2 norm of gradient
% Modification: July 30, 2020
%   penalChioce: penalty term to image, 
%       1 for L2 norm of image and 2 for L2 norm of gradient

% clear all;
% close all;
warning('off','all')
if(nargin==15)
    penalChioce = 1; % default for Tikhonov norm
end
tic;
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
rotAng = -70; % image rotation angle, (unit: degree);
flagExcludeTilt = 1; % exclude tilts
fileImgIn = [folderIn, fileName,'.tif'];
fileZernIn = [folderIn,fileName,'_coeff.txt'];
% create output folder
if isequal(exist(folderOut, 'dir'),7)
    disp(['output folder:' folderOut]);
else
    mkdir(folderOut);
    disp(['output folder created:' folderOut]);
end
disp('...Preprocessing images...');
imgs = fileIO_lvtiff2mat(fileImgIn, imgNum, repNum, cropSize, bgValue);
WriteTifStack(imgs, [folderOut 'Image_phasediversity.tif'], 32);

% % % phases and zernike coefficients
if(flagExcludeTilt==1)
    pIn = 3:zernIdxLimit; % 0: piston; 1:tilt Y; 2: tilt X;
else
    pIn = 1:zernIdxLimit;
end

% upboosting steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
% oIdx = [9 14 20 27 35 44 54];
oIdx = [10 14 20 27 35 44 54];
if(zernIdxLimit<=oIdx(1))
    zernSteps = zernIdxLimit - pIn(1) + 1;
elseif(zernIdxLimit<=oIdx(2))
    zernSteps = [oIdx(1) zernIdxLimit] - pIn(1) + 1; 
elseif(zernIdxLimit<=oIdx(3))
    zernSteps = [oIdx(1) oIdx(2) zernIdxLimit] - pIn(1) + 1; 
elseif(zernIdxLimit<=oIdx(4))
    zernSteps = [oIdx(1) oIdx(2) oIdx(3) zernIdxLimit] - pIn(1) + 1; 
elseif(zernIdxLimit<=oIdx(5))
    zernSteps = [oIdx(1) oIdx(2) oIdx(3) oIdx(4) zernIdxLimit] - pIn(1) + 1; 
else
    error('recon_zern_fileIO: Zernike coefficents order out of range');
end
% import zernike coefficients
switch(idxType)
    case 'Wyant'
        conType = 'Wyant2ANSI';
    case 'ANSI'
        conType = 'ANSI2ANSI';
    otherwise
        error('recon_zern_fileIO: wrong index type');
end
coeffsRaw = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
[phaseNum, ~] = size(coeffsRaw);
if(phaseNum<imgNum)
    error('phase number is less than image number');
end

if(flagExcludeTilt==1)
    pStart = 3; % 1:tilt Y; 2: tilt X; 3: defocus;
else
    pStart = 1;
end

pEnd = zernIdxLimit;
% unknown aberration:
coeffsInitial = coeffsRaw(1,pStart:pEnd);
% phase diversity:
coeffs_delta = coeffsRaw(2:imgNum,pStart:pEnd);
coeffs_all = [coeffsInitial;coeffs_delta;];

cTime1 = toc;
disp(['... ... time cost: ', num2str(cTime1)]);
% reconstruct zernike coefficients: estimate the unknown aberration
disp('...Reconstructing Zernike coefficients...');
[cEstimate, imgEstimate, rePar] = recon_zern(imgs, pIn, coeffs_delta, gamma,...
    alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce, rotAng);
cRMSE = rms(cEstimate(:) - coeffsInitial(:));
cTime2 = toc;
disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
save([folderOut 'data.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime2)]);

% calculate PSFs and waveFront based on input zernike coefficients
img0 = squeeze(imgs(:,:,1));
[~, PSFs, waveFronts] = gen_simu_images(img0, pIn, coeffs_all, pixelSize, ...
    lambda, NA, 'none', 10, rotAng);
WriteTifStack(PSFs, [folderOut 'PSF_phasediversity.tif'], 32);
WriteTifStack(waveFronts, [folderOut 'Wavefront_phasediversity.tif'], 32);
WriteTifStack(imgEstimate, [folderOut 'Image_estimated.tif'], 32);
[Sx, Sy] = size(img0);
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
theta = theta + deg2rad(rotAng);
waveFront = zeros(Sx, Sy, 'single');
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
WriteTifStack(waveFront, [folderOut 'Wavefront_estimated.tif'], 32);
% export zernike coefficient to txt file: imagine optic convention
switch(idxType)
    case 'Wyant'
        conType = 'ANSI2Wyant';
    case 'ANSI'
        conType = 'ANSI2ANSI';
    otherwise
        error('recon_zern_fileIO: wrong index type');
end
fileTxtOut = [folderOut, 'zernCoeffs_estimated.txt'];
coeffsOut = fileIO_coeffs2lvtxt(fileTxtOut, pIn, cEstimate, conType);

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
    savefig([folderOut 'input.fig']);
end
if(flagShowRecon == 1)
    [~, ~, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
    waveFront_original = waveFronts(:,:,1);
    wMin = min(waveFront_original(:));
    wMax = max(waveFront_original(:));
    waveFrontForShow = nan(Sx,Sx);
    waveFrontForShow(idx) = waveFront_original(idx);
    F2 = figure; subplot(2,2,1);
    pcolor(xi,xi,waveFrontForShow), shading interp
    axis square, Fc = colorbar;
    caxis([wMin wMax])
    xlabel(Fc,'\mum');
    title('Wavefront: aberrated');
    
    figure(F2); subplot(2,2,2);
    img = img0/max(img0(:));
    imshow(img,[]),colorbar;
    title('Image: aberrated');
    
    figure(F2); subplot(2,2,3);
    waveFrontForShow = nan(Sx,Sx);
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
    savefig([folderOut 'retrieval.fig']);
    
    figure; % plot Zernike coefficients
    plot(pIn,coeffsInitial,pIn,cEstimate,'LineWidth',2);
    legend( 'input','estimated');
    xlabel('Zernike Index');
    ylabel('Coeff Value');
    set(gca,'FontSize', 14);
    title('Zernike Modes');
    savefig([folderOut 'coeff.fig']);
end
cTime = toc;

