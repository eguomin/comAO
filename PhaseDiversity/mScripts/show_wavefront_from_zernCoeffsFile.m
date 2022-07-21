% show wavefront based on zernike coefficients from imagine optics
% 
% 
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;

% flieZernCoeffs = 'C:\Programs\computAO\collar_p19_coeffs_1.txt'; % Zernike coefficients file
flieZernCoeffs = 'D:\Data\20200310\collar_p19_coeffs_1.txt'; % Zernike coefficients file
coeffsRaw = importdata(flieZernCoeffs);
zernRawNum = length(coeffsRaw);
coeffsSigns = ones(1,zernRawNum); % to correct the mismatch between MATLAB and HASO
coeffsSigns([2,5, 7, 10, 12, 15]) = -1; % need further update
coeffs = coeffsRaw .* coeffsSigns;

% coeffs = coeffs * 0;
% coeffs(8) = -1;
% Imagine Opitc convetion to Fringe convention
% 1st: tilt <--> 1st: piston
pIn = 2:1+zernRawNum;
Sx = 512;
Sy = Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
% show wavefront
waveFront = nan(Sx, Sy);
waveFront(idx) = create_wavefront(pIn,coeffs,r(idx),theta(idx));
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront), shading interp
axis square, Fc = colorbar;
% caxis([-1 0.8])
xlabel(Fc,'\mum');
title('Wavefront');
