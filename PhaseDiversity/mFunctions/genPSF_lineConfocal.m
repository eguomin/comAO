function PSF = genPSF_lineConfocal(p,coeffs, Sx, pixelSize, lambda1,lambda2, NA, Sz, zStepSize, RI)
%
% By: Min Guo
% Oct. 3, 2020
Sy = Sx;
coeffs0 = zeros(size(coeffs));
PSFexc = coeffs2PSF_customAm(p,coeffs0, Sx, pixelSize, lambda1, NA, Sz,zStepSize, RI);
coeffs(p==4) = 0; % always in focus for detection path: no defocus
PSFdet = coeffs2PSF(p,coeffs, Sx, pixelSize, lambda2, NA, Sz,zStepSize, RI);

slit_width = 5; % pixels
Slit = zeros(Sx,Sy,Sz);
Slit(round(Sx/2)-floor(slit_width/2):round(Sx/2)+floor(slit_width/2),:,round(Sz/2)) = 1000;
PSF = abs(fftshift(ifftn(fftn(PSFexc).*fftn(Slit)))).*PSFdet ;