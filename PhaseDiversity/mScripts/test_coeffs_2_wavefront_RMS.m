% A simple example: convert Zernike coefficients to wavefront and calculate
% RMS parameters.

% By: Min Guo
% Last update: Aug. 12, 2021

pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;

pIn = 1:35;
% coeffs = zeros(1,32);
% coeffs(3) = 0.1;
% coeffs(5) = 0.1;
% coeffs(12) = 0.2;
% coeffs([9 12 14]) = 0.2;
coeffs = zernCoefNet_a13_0p01;
Sx = 256;
Sy = Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
% no normalization
waveFrontVector = create_wavefront(pIn,coeffs,r(idx),theta(idx));
waveFrontRMS = rms(waveFrontVector);
[waveFront0, phaseImg, staPara, puPara] = coeffs2wavefront(pIn,coeffs,Sx,...
    pixelSize, lambda, NA, 0); 

coeffs_temp = coeffs*2/staPara.rmsPhase;
[waveFront0, phaseImg, staPara, puPara] = coeffs2wavefront(pIn,coeffs_temp,Sx,...
    pixelSize, lambda, NA, 0); 
staPara.pvLength
staPara.pvPhase
staPara.rmsPhase
% show wavefront
wMin = -0.5;
wMax = 0.5;
waveFront = nan(Sx, Sy);
waveFront(idx) = waveFrontVector;

xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront');

% with normalization: RMS = root(sum(coeffs.^2))
waveFrontVector_norm = create_wavefront(pIn,coeffs,r(idx),theta(idx), 'norm');
waveFrontRMS_norm = rms(waveFrontVector_norm);
waveFrontRMS_norm2 = sqrt(sum(coeffs.^2));