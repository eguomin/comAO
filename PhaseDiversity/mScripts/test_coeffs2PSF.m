% test coeffs2PSF
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
fileFolderOut = '.\';
pIn = 1:32;
coeffs = zeros(1,32);
% % add aberrations
% coeffs(1) = -2; % tilt
% coeffs(4) = 0.5; % defocus
coeffs(5) = 0.2; % astigmatism

% % 2D PSF
Sx = 128;
PSF2D = coeffs2PSF(pIn,coeffs, Sx, pixelSize, lambda, NA);
WriteTifStack(PSF2D,[fileFolderOut 'PSF2D_astig.tif'], 32);

% % 3D PSF
Sz = 33;
zStepSize = 0.2;
RI = 1.33;
PSF3D = coeffs2PSF(pIn,coeffs, Sx, pixelSize, lambda, NA, Sz,zStepSize, RI);
WriteTifStack(PSF3D,[fileFolderOut 'PSF3D_astig.tif'], 32);