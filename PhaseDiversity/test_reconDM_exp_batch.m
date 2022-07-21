% A simple example: phase diversity algorithm on image stack (defocus diversity phase) 
% for deformable mirror calibration - batch processing.

% By: Min Guo
% Last update: July 20, 2022

% % *** default settings:
clear all;
close all;
tStart = tic;

rotAng = -90;
pixelSize = 0.08; % um
lambda = 0.532; % um
NA = 1.2;
Sx = 128;
Sy = Sx;
xi = 1: Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
theta1 = theta + deg2rad(rotAng);
% Creat the circular pupil Mask
pupilMask = zeros(Sx,Sy);
pupilMask(idx) = 1;

% create multiple defocus phase  and images
% ds = -1.6:0.4:1.6;  % defocus lenght: um
ds = -3:1:3;  % defocus lenght: um
N = length(ds);
phi_deltas = zeros(Sx, Sy, N, 'single');
% % Zernike defocus
% coeffs_delta = [];
% z = zeros(Sx, Sy);
% for i = 1:N
%     c_delta = gen_zern_coeffs(pIn, ds(i),'defocus');
%     z(idx) = create_wavefront(pIn,c_delta,r(idx),theta(idx));
%     z = length2phase*z;
%     phi_deltas(:,:,i) = z;
%     if i >=2
%         coeffs_delta = [coeffs_delta; c_delta;];
%     end
% end
% % Hanser's method
RI = 1.33;
phiunit = calc_defocusunit(Sx, pixelSize, lambda, NA, RI);
for i = 1:N
    phi_deltas(:,:,i) = ds(i)*phiunit;
end

% % reconstruction settings
flagGPU = 1;
pIn = 3:20; % 0: piston; 1:tilt X; 2: tilt Y;
zernNum = length(pIn);
itLimit = 10; % max iteration number
gamma = 1e-14; % gamma range: 1e-14 ~ 1e-6
alpha = 0;
penalChioce = 1;
% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
%     zernSteps = [14 20 27 pIn(zernNum)] - pIn(1) + 1;
zernSteps = [14 20 pIn(zernNum)] - pIn(1) + 1;
    
% % files and folders
fileFolderIn = 'D:\Data\20220718_Actuator\PD_Batch_0p03_step1um_2\processed\';
fileFolderMeasWF = 'D:\Data\20220720_Calibration_Summary\WFS_Batch_0p03\';
fileFolderPSF = [fileFolderIn 'PSF_all\'];
fileFolderWavefront = [fileFolderIn 'wavefront_all\'];
mkdir(fileFolderPSF);
mkdir(fileFolderWavefront);

fileFolderPSF_2 = [fileFolderIn 'PSF_all_-1\'];
fileFolderWavefront_2 = [fileFolderIn 'wavefront_all_-1\'];
fileFolderWavefront_net = [fileFolderIn 'wavefront_all_net\'];
mkdir(fileFolderPSF_2);
mkdir(fileFolderWavefront_2);
mkdir(fileFolderWavefront_net);
for iFile = 1:52
    disp(['image #: ' num2str(iFile)]);
    fileName = ['A_' num2str(iFile) '_1'];
    fileImgSample = [fileFolderIn fileName '.tif'];
    fileZernIn = [fileFolderMeasWF 'A_' num2str(iFile) '\A_zernCoefMeas_1.txt'];
    imgIn = single(ReadTifStack(fileImgSample));
    fileFolderOut = [fileFolderIn fileName '\'];
    fileOutSimu = [fileFolderOut 'PSF_simu_' fileName '.tif'];
    fileOutPredict = [fileFolderOut 'PSF_predict_' fileName '.tif'];
    if isequal(exist(fileFolderOut, 'dir'),7)
        disp(['output folder:' fileFolderOut]);
    else
        mkdir(fileFolderOut);
        disp(['output folder created:' fileFolderOut]);
    end
   
    % import experimental measured wavefront (by SH Sensor)
    conType = 'ANSI2ANSI';
    zernIdxLimit = 35;
    [coeffsRaw, coeffsLv] = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
    coeffs = coeffsRaw(pIn);
%     coeffs(2) = 0;  % remove defocus
    z = zeros(Sx, Sy);
    z(idx) = create_wavefront(pIn,coeffs,r(idx),theta1(idx));
    phi = 2*pi*z/lambda;
    imgsSimu = calc_multiwavefront2PSF(pupilMask, phi, phi_deltas);
    
    
    cTime1 = toc(tStart);
    % Reconstruct phase distribution
    Imgs0 = imgIn; 
    [cEstimate, ~, rePar] = recon_zern_custom(Imgs0, pIn, phi_deltas, ...
        gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);
    
    phi_c = 2*pi*cEstimate/lambda;
    % phi_c(2) = 0; % remove defocus
    phiRe = zeros(Sx, Sy);
    phiRe(idx) = create_wavefront(pIn, phi_c, r(idx), theta(idx));
    Imgs = calc_multiwavefront2PSF(pupilMask, phiRe, phi_deltas);
           
    
    phiCorrected = phi - phiRe;
    phiCorrected = angle(exp(1i*phiCorrected));
    pupilFunCorrected = pupilMask.*exp(1i*phiCorrected);
    prf = fftshift(ifft2(ifftshift(pupilFunCorrected)));
    PSF_Corrected = abs(prf).^2;
    RMS = sqrt(var(phiCorrected(idx)));
    
    WriteTifStack(imgsSimu, fileOutSimu, 32);
    WriteTifStack(Imgs, fileOutPredict, 32);
    
    save([fileFolderOut fileName '.mat']);
    cTime2 = toc(tStart);
    disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
    
    % Show results
    showWavefront = nan(Sx,Sy);
    F1 = figure; subplot(2,2,1);
    showWavefront(idx) = phi(idx);
    wMin_c = min(phi(idx));
    wMax_c = max(phi(idx));
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): SH Sensor');
    
    figure(F1), subplot(2,2,2);
    showWavefront(idx) = phiRe(idx);
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): PD estimate');
    
    figure(F1), subplot(2,2,3);
    showWavefront(idx) = phiCorrected(idx);
    pcolor(xi,xi,showWavefront), shading interp,
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): difference');
    
    savefig([fileFolderOut 'wavefront_' fileName '.fig']);
    saveas(F1,[fileFolderOut 'wavefront_' fileName '.tif'])
    saveas(F1,[fileFolderWavefront 'wavefront_' fileName '.tif'])
    
    % Generate PSF
    F2 = figure; subplot(2,2,1);
    imagesc(squeeze(imgsSimu(:,:,ceil((N+1)/2))));
    axis square, colorbar;
    title('PSF: SH simulated');
    
    figure(F2); subplot(2,2,2);
    imagesc(squeeze(Imgs(:,:,ceil((N+1)/2)))),
    axis square, colorbar;
    title('PSF: PD predicted');
    
    figure(F2); subplot(2,2,3);
    imagesc(squeeze(Imgs0(:,:,ceil((N+1)/2))));
    axis square, colorbar;
    title('raw image: focal plane');
    
    savefig([fileFolderOut 'PSF_' fileName '.fig']);
    saveas(F2,[fileFolderOut 'PSF_' fileName '.tif'])
    saveas(F2,[fileFolderPSF 'PSF_' fileName '.tif'])
    close all;
    
    coeffs_1 = coeffs;
    phi_1 = phi;
    phi_c_1 = phi_c;
    phiRe_1 = phiRe;
    
    
    fileName = ['A_' num2str(iFile) '_-1'];
    fileImgSample = [fileFolderIn fileName '.tif'];
    fileZernIn = [fileFolderMeasWF 'A_' num2str(iFile) '\A_zernCoefMeas_-1.txt'];
    imgIn = single(ReadTifStack(fileImgSample));
    fileFolderOut = [fileFolderIn fileName '\'];
    fileOutSimu = [fileFolderOut 'PSF_simu_' fileName '.tif'];
    fileOutPredict = [fileFolderOut 'PSF_predict_' fileName '.tif'];
    if isequal(exist(fileFolderOut, 'dir'),7)
        disp(['output folder:' fileFolderOut]);
    else
        mkdir(fileFolderOut);
        disp(['output folder created:' fileFolderOut]);
    end
   
    % import experimental measured wavefront (by SH Sensor)
    conType = 'ANSI2ANSI';
    zernIdxLimit = 35;
    [coeffsRaw, coeffsLv] = fileIO_lvtxt2idx(fileZernIn, conType, zernIdxLimit);
    coeffs = coeffsRaw(pIn);
%     coeffs(2) = 0;  % remove defocus
    z = zeros(Sx, Sy);
    z(idx) = create_wavefront(pIn,coeffs,r(idx),theta1(idx));
    phi = 2*pi*z/lambda;
    imgsSimu = calc_multiwavefront2PSF(pupilMask, phi, phi_deltas);
    
    cTime1 = toc(tStart);
    % Reconstruct phase distribution
    Imgs0 = imgIn; 
    [cEstimate, ~, rePar] = recon_zern_custom(Imgs0, pIn, phi_deltas, ...
        gamma, alpha, itLimit, zernSteps, pixelSize, lambda, NA, flagGPU, penalChioce);

     phi_c = 2*pi*cEstimate/lambda;
    % phi_c(2) = 0; % remove defocus
    phiRe = zeros(Sx, Sy);
    phiRe(idx) = create_wavefront(pIn, phi_c, r(idx), theta(idx));
    Imgs = calc_multiwavefront2PSF(pupilMask, phiRe, phi_deltas);
           
    
    phiCorrected = phi - phiRe;
    phiCorrected = angle(exp(1i*phiCorrected));
    pupilFunCorrected = pupilMask.*exp(1i*phiCorrected);
    prf = fftshift(ifft2(ifftshift(pupilFunCorrected)));
    PSF_Corrected = abs(prf).^2;
    RMS = sqrt(var(phiCorrected(idx)));
    
    WriteTifStack(imgsSimu, fileOutSimu, 32);
    WriteTifStack(Imgs, fileOutPredict, 32);
    
    save([fileFolderOut fileName '.mat']);
    cTime2 = toc(tStart);
    disp(['... ... time cost: ', num2str(cTime2-cTime1)]);
    
    % Show results
    showWavefront = nan(Sx,Sy);
    F1 = figure; subplot(2,2,1);
    showWavefront(idx) = phi(idx);
    wMin_c = min(phi(idx));
    wMax_c = max(phi(idx));
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): SH Sensor');
    
    figure(F1), subplot(2,2,2);
    showWavefront(idx) = phiRe(idx);
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): PD estimate');
    
    figure(F1), subplot(2,2,3);
    showWavefront(idx) = phiCorrected(idx);
    pcolor(xi,xi,showWavefront), shading interp,
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): difference');
    
    savefig([fileFolderOut 'wavefront_' fileName '.fig']);
    saveas(F1,[fileFolderOut 'wavefront_' fileName '.tif'])
    saveas(F1,[fileFolderWavefront_2 'wavefront_' fileName '.tif'])
    
    % Generate PSF
    F2 = figure; subplot(2,2,1);
    imagesc(squeeze(imgsSimu(:,:,ceil((N+1)/2))));
    axis square, colorbar;
    title('PSF: SH simulated');
    
    figure(F2); subplot(2,2,2);
    imagesc(squeeze(Imgs(:,:,ceil((N+1)/2)))),
    axis square, colorbar;
    title('PSF: PD predicted');
    
    figure(F2); subplot(2,2,3);
    imagesc(squeeze(Imgs0(:,:,ceil((N+1)/2))));
    axis square, colorbar;
    title('raw image: focal plane');
    
    savefig([fileFolderOut 'PSF_' fileName '.fig']);
    saveas(F2,[fileFolderOut 'PSF_' fileName '.tif'])
    saveas(F2,[fileFolderPSF_2 'PSF_' fileName '.tif'])
    close all;
    
    coeffs_2 = coeffs;
    phi_2 = phi;
    phi_c_2 = phi_c;
    phiRe_2 = phiRe;
    
    coeffs = coeffs_1 - coeffs_2;
    phi = phi_1 - phi_2;
    phi_c = phi_c_1 - phi_c_2;
    phiRe = phiRe_1 - phiRe_2;
    
    % Show results
    showWavefront = nan(Sx,Sy);
    F1 = figure; subplot(2,2,1);
    showWavefront(idx) = phi(idx);
    wMin_c = min(phi(idx));
    wMax_c = max(phi(idx));
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): SH Sensor');
    
    figure(F1), subplot(2,2,2);
    showWavefront(idx) = phiRe(idx);
    pcolor(xi,xi,showWavefront), shading interp
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): PD estimate');
    
    figure(F1), subplot(2,2,3);
    phiCorrected = phi - phiRe;
    showWavefront(idx) = phiCorrected(idx);
    pcolor(xi,xi,showWavefront), shading interp,
    axis square, Fc = colorbar;
    caxis([wMin_c wMax_c])
    xlabel(Fc,'rad');
    title('{\phi}(r,\theta): difference');
    
    savefig([fileFolderOut 'wavefront_a_' num2str(iFile) '.fig']);
    saveas(F1,[fileFolderOut 'wavefront_a_' num2str(iFile) '.tif'])
    saveas(F1,[fileFolderWavefront_net 'wavefront_a_' num2str(iFile) '.tif'])
    
    close all;
end
