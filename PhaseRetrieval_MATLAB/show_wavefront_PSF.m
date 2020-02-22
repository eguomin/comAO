% show wavefront and PSF
pixelSize = 0.096; % um
lambda = 0.550; % um
NA = 1.2;

% pIn = 1:32;
% cEstimate = zeros(1,32);
Sx = 768;
Sy = Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
% show wavefront
wMin = -0.2;
wMax = 0.2;
waveFront = nan(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
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
z(idx) = create_wavefront(pIn,cEstimate,r(idx),theta(idx));
phi = z*length2phase;
pupilFun = pupilMask.*exp(1i*phi);
prf = fftshift(ifft2(ifftshift(pupilFun)));
PSF = abs(prf).^2;
figure;
imshow(PSF/max(PSF(:)),[]),colorbar;
title('PSF');
WriteTifStack(PSF, 'PSF_calculated.tif', 32);