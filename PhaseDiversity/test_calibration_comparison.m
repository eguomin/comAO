% show calibration comparison: PD and WFS

% Min Guo
% July 20, 2022

clear all
close all;
pIn = 3:20; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);

pixelSize = 0.2; % um, larger pixel size to make pupil larger in display
lambda = 0.532; % um
NA = 1.2;
Sx = 128;
Sy = Sx;
xi = 1: Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);

% path interaction matrix (influence matrix):
fileFolderMain = 'D:\Data\20220720_Calibration_Summary\Calibration_Check\';
fileFolderMeasWF = [fileFolderMain 'zernCheck_WFS_0p3\'];
fileFolderPDWF = [fileFolderMain 'zernCheck_PD_0p3\'];
fileFolderOut = [fileFolderMain 'analysis\'];
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

fileFolderOutInterWF = [fileFolderOut 'interWF\'];
if isequal(exist(fileFolderOutInterWF, 'dir'),7)
    disp(['output folder:' fileFolderOutInterWF]);
else
    mkdir(fileFolderOutInterWF);
    disp(['output folder created:' fileFolderOutInterWF]);
end

zernAm = 0.3; %um
% import interaction matrix: WFS [zernNum x actNum]
zernMatrixGT = zeros(zernNum, zernNum);
zernMatrixWFS = zeros(zernNum, zernNum);
zernMatrixPD = zeros(zernNum, zernNum);
conType = 'ANSI2ANSI';
showWavefront = nan(Sx,Sy);
for i = 1: zernNum
    iFile = pIn(i);
    disp(['zern #: ' num2str(iFile)]);
    coeffs = zeros(1, zernNum);
    coeffs(i) = zernAm;
    zernMatrixGT(:,i) = coeffs';
    z = create_wavefront(pIn,coeffs,r(idx),theta(idx));
    wMin_c = min(z);
    wMax_c = max(z);
    F1 = figure; subplot(2,2,1);
    showWavefront(idx) = z;
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'\mum');
    title('{\phi}(r,\theta): Ground truth');
    
    fileZernIn = [fileFolderMeasWF 'zernCheck_WFS_0p3_Z_' num2str(iFile) '\Z_zernCoefNet.txt'];
    zernIdxLimit = 35;
    [coeffsRaw, ~] = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
    coeffs = coeffsRaw(pIn);
    zernMatrixWFS(:,i) = coeffs';
    z = create_wavefront(pIn,coeffs,r(idx),theta(idx));
    figure(F1), subplot(2,2,3);
    showWavefront(idx) = z;
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'\mum');
    title('{\phi}(r,\theta): SH Sensor');
    
    fileZernIn = [fileFolderPDWF 'zernCheck_PD_0p3_Z_' num2str(iFile) '\Z_zernCoefNet.txt'];
    zernIdxLimit = 35;
    [coeffsRaw, coeffsLv] = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
    coeffs = coeffsRaw(pIn);
    zernMatrixPD(:,i) = coeffs';
    z = create_wavefront(pIn,coeffs,r(idx),theta(idx));
    figure(F1), subplot(2,2,4);
    showWavefront(idx) = z;
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'\mum');
    title('{\phi}(r,\theta): PD estimate');
    saveas(F1,[fileFolderOutInterWF 'wavefront_a_' num2str(i) '.tif'])
    close all;
end
save([fileFolderOut 'calibration_results.mat']);

CLIM = [-0.5 zernAm];
F2 = figure;
imagesc(flipud(zernMatrixGT),CLIM);
axis square, Fc = colorbar;
xlabel(Fc,'Amplitude');
saveas(F2,[fileFolderOut 'ZernikeMap_GT.tif'])

F3 = figure;
imagesc(flipud(zernMatrixWFS),CLIM);
axis square, Fc = colorbar;
xlabel(Fc,'Amplitude');
saveas(F3,[fileFolderOut 'ZernikeMap_WFS.tif'])

F4 = figure;
imagesc(flipud(zernMatrixPD),CLIM);
axis square, Fc = colorbar;
xlabel(Fc,'Amplitude');
saveas(F4,[fileFolderOut 'ZernikeMap_PD.tif'])
