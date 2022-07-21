% A simple example: convert Zernike coefficients to 2P PSF (camera detection),
%   and show the images.

% By: Min Guo
% Last update: Oct. 21, 2020
close all
pixelSize = 0.096; % um
lambdaDet = 0.532; % um
lambdaExc = 0.800; % um
NA = 1.2;

normFlag = 0;
% flagPupilMatch = 1;
fileFolderOut = 'C:\Programs\computAO\Data\PSF2P\astig_0p2\';
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
% abType1 = 'astig';
abType1 = 'astig';
abRMS1 = 0.2;
abValue1 = abRMS1 * sqrt(5); % defocus: Z3 = sqrt(3), astig: Z5 = sqrt(6)
            % coma: Z8 = sqrt(8), trefoil: Z9 = sqrt(8)
            % sphe: Z12 = sqrt(5)
            % idx: customize aberrations
coeffs = gen_zern_coeffs(pIn,abValue1,abType1);
% coeffs([9 12 14]) = 0.2;
Sx = 256;
Sy = Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambdaDet, NA);
waveFrontVector = create_wavefront(pIn,coeffs,r(idx),theta(idx));
% show detection wavefront
wMin = -0.3;
wMax = 0.3;
waveFront = nan(Sx, Sy);
waveFront(idx) = waveFrontVector;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Det Wavefront');

[r2, theta2, idx2] = def_pupilcoor(Sx, pixelSize, lambdaExc, NA);
waveFront2 = nan(Sx, Sy);
waveFront2(idx2) = waveFront(idx2);
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront2), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Excitation Wavefront');

[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambdaExc, NA);
waveFrontVector3 = create_wavefront(pIn,coeffs,r(idx),theta(idx));
% show detection wavefront
wMin = -0.3;
wMax = 0.3;
waveFront3 = nan(Sx, Sy);
waveFront3(idx) = waveFrontVector3;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront3), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Excitation Wavefront');

% detection PSF
% % 3D PSF
Sz = Sx;
zStepSize = pixelSize;
RI = 1.33;
Sox = round((Sx+1)/2);
Soy = round((Sy+1)/2);
Soz = round((Sz+1)/2);

flagPupilMatch = 0;
[PSF_2P, PSF_Det, PSF_Exc, PSF_Exc2P_ave] = genPSF_Camera2P(pIn,coeffs, Sx, ...
    pixelSize, lambdaExc,lambdaDet, NA, Sz, zStepSize, RI, normFlag, flagPupilMatch);
PSF_2P_line0 = PSF_2P(Sox,Soy,:);
PSF_Exc2P_ave_line0 = PSF_Exc2P_ave(Sox,Soy,:);

flagPupilMatch = 1;
[PSF_2P, PSF_Det, PSF_Exc, PSF_Exc2P_ave] = genPSF_Camera2P(pIn,coeffs, Sx, ...
    pixelSize, lambdaExc,lambdaDet, NA, Sz, zStepSize, RI, normFlag, flagPupilMatch);
PSF_2P_line1 = PSF_2P(Sox,Soy,:);
PSF_Exc2P_ave_line1 = PSF_Exc2P_ave(Sox,Soy,:);

xi =1:Sz;
xi = (xi- Soz)*pixelSize;
PSF_Exc2P_ave_line0 = PSF_Exc2P_ave_line0(:)/max(PSF_Exc2P_ave_line0(:));
PSF_Exc2P_ave_line1 = PSF_Exc2P_ave_line1(:)/max(PSF_Exc2P_ave_line1(:));
figure, plot(xi, PSF_Exc2P_ave_line0, ...
    xi, PSF_Exc2P_ave_line1, 'LineWidth',2);
ylim([0 1.1]);
xlim([-5 5]);
legend('Pupil model 1', 'Pupil model 2');
xlabel('z positions ({\mu}m)');
ylabel('Normalized Intensity');
title('Z Profiles of Exc2P Ave PSF');

PSF_2P_line0 = PSF_2P_line0(:)/max(PSF_2P_line0(:));
PSF_2P_line1 = PSF_2P_line1(:)/max(PSF_2P_line1(:));
figure, plot(xi, PSF_2P_line0, ...
    xi, PSF_2P_line1, 'LineWidth',2);
ylim([0 1.1]);
xlim([-5 5]);
legend('Pupil model 1', 'Pupil model 2');
xlabel('z positions ({\mu}m)');
ylabel('Normalized Intensity');
title('Z Profiles of final 2P PSF');

