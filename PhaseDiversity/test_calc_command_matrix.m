% calculate command matrix, and then generate commands based on this command matrix

% Min Guo
% July 20, 2022

clear all
close all;
actNum = 52;
pIn = 3:20; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);

% path interaction matrix (influence matrix):
fileFolderMeasWF = 'D:\Data\20220720_Calibration_Summary\WFS_Batch_0p03\';
fileFolderPDWF = 'D:\Data\20220720_Calibration_Summary\PD_Batch_0p03_step0p5um\processed\';
fileFolderOut = 'D:\Data\20220720_Calibration_Summary\Calc_PD_Batch_0p03_step0p5\';
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

volAm = 0.03;
% import interaction matrix: WFS [zernNum x actNum]
interMatrixWFS = zeros(zernNum, actNum);
conType = 'ANSI2ANSI';
for iFile = 1:actNum
    disp(['image #: ' num2str(iFile)]);
    fileZernIn = [fileFolderMeasWF 'A_' num2str(iFile) '\A_zernCoefMeas_1.txt'];
    zernIdxLimit = 35;
    [coeffsRaw, coeffsLv] = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
    coeffs_1 = coeffsRaw(pIn);
    fileZernIn = [fileFolderMeasWF 'A_' num2str(iFile) '\A_zernCoefMeas_-1.txt'];
    zernIdxLimit = 35;
    [coeffsRaw, coeffsLv] = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
    coeffs_2 = coeffsRaw(pIn);
    coeffs = (coeffs_1 - coeffs_2)/volAm/2;
    coeffs = coeffs';
    interMatrixWFS(:, iFile) = coeffs;
end

% import interaction matrix: PD [zernNum x actNum]
interMatrixPD = zeros(zernNum, actNum);
rotAng = 90;
lambda = 0.532; % um
for iFile = 1:actNum
    disp(['image #: ' num2str(iFile)]);
    fileName = ['A_' num2str(iFile) '_1'];
    fileZernIn = [fileFolderPDWF fileName '\' fileName '.mat'];
    load(fileZernIn, 'phi_c') 
    coeffs_1 = phi_c * lambda/pi/2;
    fileName = ['A_' num2str(iFile) '_-1'];
    fileZernIn = [fileFolderPDWF fileName '\' fileName '.mat'];
    load(fileZernIn, 'phi_c') 
    coeffs_2 = phi_c * lambda/pi/2;
    coeffs = (coeffs_1 - coeffs_2)/volAm/2;
    % rotate wavefront to match WFS measurement
    coeffs = coeffs_rot(coeffs, pIn, rotAng);
    coeffs = coeffs';
    interMatrixPD(:, iFile) = coeffs;
end

% compute command matrix: [actNum x zernNum]
threValue = 0.005; % from Martin paper
cmdMatrixWFS = pinv(interMatrixWFS, threValue);
cmdMatrixPD = pinv(interMatrixPD, threValue);

% compute command for each Zernike component: 0.3 um amplitude
zernAm = 0.3; % um
zernCmdWFS = cmdMatrixWFS * zernAm;
zernCmdPD = cmdMatrixPD * zernAm;

% % compute command for diversity phases: 4 phases consistent with LabVIEW
% settings
zernDivCoeffs = zeros(zernNum, 4);
zernDivCoeffs(1,1) = 0.5;
zernDivCoeffs(2,1) = 0.2;
zernDivCoeffs(2,2) = -0.2;
zernDivCoeffs(3,2) = 0.5;
zernDivCoeffs(1,3) = 0.5;
zernDivCoeffs(2,3) = -0.2;
zernDivCoeffs(2,4) = 0.2;
zernDivCoeffs(3,4) = 0.5;
zernDivCmdWFS = cmdMatrixWFS * zernDivCoeffs;
zernDivCmdPD = cmdMatrixPD * zernDivCoeffs;

% save files
save([fileFolderOut 'calc_results.mat']);
dlmwrite([fileFolderOut 'interMatrixWFS.dat'], interMatrixWFS ,'delimiter','\t','precision','%.6f');
dlmwrite([fileFolderOut 'interMatrixPD.dat'], interMatrixPD ,'delimiter','\t','precision','%.6f');
dlmwrite([fileFolderOut 'cmdMatrixWFS.dat'], cmdMatrixWFS ,'delimiter','\t','precision','%.6f');
dlmwrite([fileFolderOut 'cmdMatrixPD.dat'], cmdMatrixPD ,'delimiter','\t','precision','%.6f');
for i = 1:zernNum
    dlmwrite([fileFolderOut 'cmd_WFS_zern_' num2str(pIn(i)) '.dat'], ...
        zernCmdWFS(:,i)' ,'delimiter','\t','precision','%.6f');
    dlmwrite([fileFolderOut 'cmd_PD_zern_' num2str(pIn(i)) '.dat'], ...
        zernCmdPD(:,i)' ,'delimiter','\t','precision','%.6f');
end
dlmwrite([fileFolderOut 'zernDivCmdWFS.dat'], zernDivCmdWFS','delimiter','\t','precision','%.6f');
dlmwrite([fileFolderOut 'zernDivCmdPD.dat'], zernDivCmdPD','delimiter','\t','precision','%.6f');

% show wavefront
pixelSize = 0.2; % um, larger pixel size to make pupil larger in display
lambda = 0.532; % um
NA = 1.2;
Sx = 128;
Sy = Sx;
xi = 1: Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
fileFolderOutInterWF = [fileFolderOut 'interWF\'];
if isequal(exist(fileFolderOutInterWF, 'dir'),7)
    disp(['output folder:' fileFolderOutInterWF]);
else
    mkdir(fileFolderOutInterWF);
    disp(['output folder created:' fileFolderOutInterWF]);
end
lambda = 0.532; % um
for i = 1:actNum
    coeffs = interMatrixWFS(:,i)';
    coeffs = coeffs_rot(coeffs, pIn, -rotAng); % rotate just for display
    z = create_wavefront(pIn,coeffs,r(idx),theta(idx));
    % z = z*2*pi/lambda;
    showWavefront = nan(Sx,Sy);
    F1 = figure; subplot(2,2,1);
    showWavefront(idx) = z;
    wMin_c = min(z);
    wMax_c = max(z);
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'\mum/V');
    title('{\phi}(r,\theta): SH Sensor');
    z1 = z;
    
    coeffs = interMatrixPD(:,i)';
    coeffs = coeffs_rot(coeffs, pIn, -rotAng);
    z = create_wavefront(pIn,coeffs,r(idx),theta(idx));
    % z = z*2*pi/lambda;
    figure(F1), subplot(2,2,2);
    showWavefront(idx) = z;
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'\mum/V');
    title('{\phi}(r,\theta): PD estimate');
    z2 = z;
    
    z = z1-z2;
    figure(F1), subplot(2,2,4);
    showWavefront(idx) = z;
    pcolor(xi,xi,showWavefront), shading interp,
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'\mum/V');
    title('{\phi}(r,\theta): difference');
    saveas(F1,[fileFolderOutInterWF 'wavefront_a_' num2str(i) '.tif'])
    close all;
end

