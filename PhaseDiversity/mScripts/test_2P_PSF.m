% A simple example: convert Zernike coefficients to 2P PSF (camera detection),
%   and show the images.

% By: Min Guo
% Last update: Oct. 21, 2020
% close all
pixelSize = 0.096; % um
lambdaDet = 0.532; % um
lambdaExc = 0.800; % um
NA = 1.2;

normFlag = 0;
flagPupilMatch = 0;
fileFolderOut = 'C:\Programs\computAO\Data\PSF2P\coma_0p2_unmatchPupil\';
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
% abType1 = 'astig';
abType1 = 'coma';
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
waveFrontRMS = rms(waveFrontVector);
% show detection wavefront
wMin = -0.3;
wMax = 0.3;
waveFront = nan(Sx, Sy);
waveFront(idx) = waveFrontVector;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront2), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront');

% detection PSF
% % 3D PSF
Sz = Sx;
zStepSize = pixelSize;
RI = 1.33;


[PSF_2P, PSF_Det, PSF_Exc, PSF_Exc2P_ave] = genPSF_Camera2P(pIn,coeffs, Sx, ...
    pixelSize, lambdaExc,lambdaDet, NA, Sz, zStepSize, RI, normFlag, flagPupilMatch);
Sox = round((Sx+1)/2);
Soy = round((Sy+1)/2);
Soz = round((Sz+1)/2);
PSF_Det_line = PSF_Det(Sox,Soy,:);
PSF_Det_line = PSF_Det_line(:)/max(PSF_Det_line(:));
PSF_2P_line = PSF_2P(Sox,Soy,:);
PSF_2P_line = PSF_2P_line(:)/max(PSF_2P_line(:));
PSF_Exc2P_ave_line = PSF_Exc2P_ave(Sox,Soy,:);
PSF_Exc2P_ave_line = PSF_Exc2P_ave_line(:)/max(PSF_Exc2P_ave_line(:));
xi =1:Sz;
xi = (xi- Soz)*pixelSize;
figure, plot(xi, PSF_Det_line, xi, PSF_Exc2P_ave_line, ...
    xi, PSF_2P_line, 'LineWidth',2);
ylim([0 1.1]);
xlim([-12 12]);
legend('Det-widefield', 'Exc2P Ave', 'final 2P PSF');
title('Z Profiles of PSFs');
WriteTifStack(PSF_Det,[fileFolderOut 'PSF3D_Det.tif'], 32);
WriteTifStack(PSF_Exc,[fileFolderOut, 'PSF3D_Exc.tif'], 32);
WriteTifStack(PSF_Exc2P_ave, [fileFolderOut, 'PSF3D_Exc2P_ave.tif'], 32);
WriteTifStack(PSF_2P, [fileFolderOut, 'PSF3D_2P.tif'], 32);

