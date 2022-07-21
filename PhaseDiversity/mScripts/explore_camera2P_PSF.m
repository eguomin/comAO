% explore 2P PSF 
% A simple example: convert Zernike coefficients to 2P PSF (camera detection),
%   and show the images.

% By: Min Guo
% Last update: Oct. 21, 2020
% close all
pixelSize = 0.096; % um
lambdaDet = 0.532; % um
lambdaExc = 0.800; % um
NA = 1.2;
Sx = 257;
Sy = Sx;
Sz = Sx;
zStepSize = pixelSize;
RI = 1.33;
normFlag = 0; % 1: normalize maximum to 1; 0: no normalization;
flagPupilMatch = 1; % 1: match pupil size of excitation wavelength; 0: unmatch
fileFolderOut = 'C:\Programs\computAO\Data\PSF2P\explore_sphe\';
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
% abType1 = 'astig';
abType1 = 'sphe';
abRMS1 = 1;
abValue1 = abRMS1 * sqrt(5); % defocus: Z3 = sqrt(3), astig: Z5 = sqrt(6)
            % coma: Z8 = sqrt(8), trefoil: Z9 = sqrt(8)
            % sphe: Z12 = sqrt(5)
            % idx: customize aberrations
coeffs = gen_zern_coeffs(pIn,abValue1,abType1);
% coeffs([9 12 14]) = 0.2;

amRatios = 0:0.1:1.2;
rNum = length(amRatios);
Sox = round((Sx+1)/2);
Soy = round((Sy+1)/2);
Soz = round((Sz+1)/2);
PSF_2P_lines = zeros(rNum,Sz);
PSF_Det_lines = zeros(rNum,Sz);
PSF_Exc2P_ave_lines = zeros(rNum,Sz);
psfSizePSF_2P = zeros(rNum,3);
psfSizePSF_Det = zeros(rNum,3);
psfSizePSF_Exc = zeros(rNum,3);
psfSizePSF_Exc2P_ave = zeros(rNum,1);
psfMaxPSF_2P = zeros(rNum,1);
psfMaxPSF_Det = zeros(rNum,1);
psfMaxPSF_Exc2P_ave = zeros(rNum,1);
psfMaxPSF_2P_central = zeros(rNum,1);
psfMaxPSF_Det_central = zeros(rNum,1);
for j = 1:rNum 
    coeffs0 = coeffs * amRatios(j);
    [PSF_2P, PSF_Det, PSF_Exc, PSF_Exc2P_ave] = genPSF_camera2P(pIn,coeffs0, Sx, ...
        pixelSize, lambdaExc,lambdaDet, NA, Sz, zStepSize, RI, normFlag, flagPupilMatch);
    PSF_2P_line = PSF_2P(Sox,Soy,:);
    PSF_2P_lines(j,:) = PSF_2P_line(:)/max(PSF_2P_line(:));
    PSF_Det_line = PSF_Det(Sox,Soy,:);
    PSF_Det_lines(j,:) = PSF_Det_line(:)/max(PSF_Det_line(:));
    PSF_Exc2P_ave_line = PSF_Exc2P_ave(Sox,Soy,:);
    PSF_Exc2P_ave_lines(j,:) = PSF_Exc2P_ave_line(:)/max(PSF_Exc2P_ave_line(:));
%     psfSizePSF_2P(j,:) = fwhm_PSF(PSF_2P,pixelSize);
%     psfSizePSF_Det(j,:) = fwhm_PSF(PSF_Det,pixelSize);
%     psfSizePSF_Exc(j,:) = fwhm_PSF(PSF_Exc,pixelSize);
%     psfSizePSF_Exc2P_ave(j) = fwhm(1:Sz,PSF_Exc2P_ave_line(:))*pixelSize;
    psfMaxPSF_2P(j) = max(PSF_2P(:));
    psfMaxPSF_Det(j) = max(PSF_Det(:));
    psfMaxPSF_Exc2P_ave(j) = max(PSF_Exc2P_ave_line(:));
    centralPlane = PSF_2P(:,:,Soz);
    psfMaxPSF_2P_central(j) = max(centralPlane(:));
    centralPlane = PSF_Det(:,:,Soz);
    psfMaxPSF_Det_central(j) = max(centralPlane(:));
end

psfMaxPSF_2P_norm = psfMaxPSF_2P/max(psfMaxPSF_2P(:));
psfMaxPSF_Det_norm = psfMaxPSF_Det/max(psfMaxPSF_Det(:));
psfMaxPSF_Exc2P_ave_norm = psfMaxPSF_Exc2P_ave/max(psfMaxPSF_Exc2P_ave(:));
figure, plot(amRatios, psfMaxPSF_Exc2P_ave_norm, amRatios, psfMaxPSF_Det_norm, ...
 amRatios, psfMaxPSF_2P_norm,'LineWidth',2);
ylim([0 1.1]);
legend('Exc2P ave', 'Det-widefield', 'final 2P');
xlabel('Aberration RMS ({\mu}m)');
ylabel('Normalized Peak Intensity');
title('Mamixum intensity (normalized to aberration-free case)');

psfMaxPSF_2P_central_norm = psfMaxPSF_2P_central/max(psfMaxPSF_2P_central(:));
psfMaxPSF_Det_central_norm = psfMaxPSF_Det_central/max(psfMaxPSF_Det_central(:));
figure, plot(amRatios, psfMaxPSF_Det_central_norm, ...
 amRatios, psfMaxPSF_2P_central_norm,'LineWidth',2);
ylim([0 1.1]);
legend('Det-widefield', 'final 2P');
xlabel('Aberration RMS ({\mu}m)');
ylabel('Normalized Peak Intensity');
title('Mamixum intensity at focal plane(normalized to aberration-free case)');

figure, plot(amRatios, psfSizePSF_Exc2P_ave, amRatios, psfSizePSF_Det(:,3),...
 amRatios, psfSizePSF_2P(:,3), 'LineWidth',2);
legend('Exc2P ave', 'Det-widefield', 'final 2P');
xlabel('Aberration RMS ({\mu}m)');
ylabel('FWHM z ({\mu}m)');
title('FWHM of PSFs (z direction)');

xi =1:Sz;
xi = (xi- Soz)*pixelSize;
figure, plot(xi, PSF_Exc2P_ave_lines(1,:), xi, PSF_Exc2P_ave_lines(3,:), ...
    xi, PSF_Exc2P_ave_lines(5,:), xi, PSF_Exc2P_ave_lines(rNum,:), 'LineWidth',2);
ylim([0 1.1]);
xlim([-12 12]);
legend(['RMS=', num2str(amRatios(1))], ['RMS=', num2str(amRatios(3))], ...
['RMS=', num2str(amRatios(5))], ['RMS=', num2str(amRatios(rNum))]);
xlabel('z positions ({\mu}m)');
ylabel('Exc2P ave intensity (normalized)');
title('Z Profiles of average of 2P exc PSF');

i = 1;
figure, plot(xi, PSF_Exc2P_ave_lines(i,:), xi, PSF_Det_lines(i,:), ...
    xi, PSF_2P_lines(i,:), 'LineWidth',2);
ylim([0 1.1]);
xlim([-12 12]);
legend('Exc2P ave', 'Det-widefield', 'final 2P');
xlabel('z positions ({\mu}m)');
ylabel('Normalized Intensity');
title('Z Profiles of PSFs: Aberration - free');

i = 3;
figure, plot(xi, PSF_Exc2P_ave_lines(i,:), xi, PSF_Det_lines(i,:), ...
    xi, PSF_2P_lines(i,:), 'LineWidth',2);
ylim([0 1.1]);
xlim([-12 12]);
legend('Exc2P ave', 'Det-widefield', 'final 2P');
xlabel('z positions ({\mu}m)');
ylabel('Normalized Intensity');
title(['Z Profiles of PSFs: RMS=', num2str(amRatios(i))]);

i = 5;
figure, plot(xi, PSF_Exc2P_ave_lines(i,:), xi, PSF_Det_lines(i,:), ...
    xi, PSF_2P_lines(i,:), 'LineWidth',2);
ylim([0 1.1]);
xlim([-12 12]);
legend('Exc2P ave', 'Det-widefield', 'final 2P');
xlabel('z positions ({\mu}m)');
ylabel('Normalized Intensity');
title(['Z Profiles of PSFs: RMS=', num2str(amRatios(i))]);

save([fileFolderOut, 'matlabData.mat']); 

