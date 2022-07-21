close all
clear all

fileFolderIn = 'D:\Data\20220415_Actuator\bead_0p5um_2\processed\';
fileFolderOut = [fileFolderIn 'results_net\'];
fileFolderWavefront = [fileFolderOut 'wavefront_all\'];
mkdir(fileFolderOut);
mkdir(fileFolderWavefront);
fileName = 'a_7_0p05';
fileName2 = 'a_7_-0p05';
load([fileFolderIn fileName '\data.mat'], 'coeffs', 'phi', 'phi_c', 'phiRe', 'Sx', 'Sy', 'idx', 'xi');
coeffs_1 = coeffs;
phi_1 = phi;
phi_c_1 = phi_c;
phiRe_1 = phiRe;

load([fileFolderIn fileName2 '\data.mat'], 'coeffs', 'phi', 'phi_c', 'phiRe');
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
F1 = figure; subplot(1,3,1);
showWavefront(idx) = phi(idx);
wMin_c = min(phi(idx));
wMax_c = max(phi(idx));
pcolor(xi,xi,showWavefront), shading interp
axis square, Fc = colorbar;
caxis([wMin_c wMax_c])
xlabel(Fc,'rad');
title('{\phi}(r,\theta): SH Sensor');

figure(F1), subplot(1,3,2);
showWavefront(idx) = phiRe(idx);
pcolor(xi,xi,showWavefront), shading interp
axis square, Fc = colorbar;
% caxis([wMin_c wMax_c])
xlabel(Fc,'rad');
title('{\phi}(r,\theta): PD estimate');

figure(F1), subplot(1,3,3);
phiCorrected = phi - phiRe;
showWavefront(idx) = phiCorrected(idx);
pcolor(xi,xi,showWavefront), shading interp,
axis square, Fc = colorbar;
caxis([wMin_c wMax_c])
xlabel(Fc,'rad');
title('{\phi}(r,\theta): difference');

savefig([fileFolderOut 'wavefront_' fileName '.fig']);
saveas(F1,[fileFolderOut 'wavefront_' fileName '.tif'])
saveas(F1,[fileFolderWavefront 'wavefront_' fileName '.tif'])
