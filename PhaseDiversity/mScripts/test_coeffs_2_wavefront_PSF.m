% A simple example: convert Zernike coefficients to wavefront and PSF,
%   and show the images.

% By: Min Guo
% Last update: Sep. 1, 2020

pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;

pIn = 1:32;
coeffs = zeros(1,32);
coeffs(3) = 0.3;
coeffs(6) = -0.3;
coeffs(8) = 0.5;
% coeffs([9 12 14]) = 0.2;
Sx = 512;
Sy = Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFrontVector = create_wavefront(pIn,coeffs,r(idx),theta(idx));
waveFrontRMS = rms(waveFrontVector);
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

% get PSF
pupilMask = zeros(Sx, Sy);
pupilMask(idx) = 1;
length2phase = 2*pi/lambda;
z = zeros(Sx, Sy);
z(idx) = waveFrontVector;
phi = z*length2phase;
pupilFun = pupilMask.*exp(1i*phi);
prf = fftshift(ifft2(ifftshift(pupilFun)));
PSF = abs(prf).^2;
figure;
Sa = 40;
PSF2 = PSF(round((Sx+1)/2)-Sa:round((Sx+1)/2)+Sa,round((Sy+1)/2)-Sa:round((Sy+1)/2)+Sa);
imshow(PSF2/max(PSF2(:)),[]),colorbar;
title('PSF');