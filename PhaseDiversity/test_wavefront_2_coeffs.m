% A simple example: convert Zernike coefficients to wavefront, and convert
% back to Zernike coefficients

% By: Min Guo
% Last update: July 19, 2022

pixelSize = 0.080; % um
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

coeffsOut = wavefront2coeffs(waveFrontVector , pIn, r(idx),theta(idx));
waveFrontVector = create_wavefront(pIn,coeffsOut,r(idx),theta(idx));
% show wavefront
waveFrontOut = nan(Sx, Sy);
waveFrontOut(idx) = waveFrontVector;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFrontOut), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Wavefront');