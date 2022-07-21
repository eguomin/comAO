% show wavefront based on zernike coefficients from imagine optics
% 
% 
close all
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
zNum = 20;
flieZernCoeffs1 = ['C:\Programs\computAO\Data\AbeComponentCheck\20200711\ZernIdx_', ...
    num2str(zNum), '_diff.txt']; % Zernike coefficients file
flieZernCoeffs2 = ['C:\Programs\computAO\Data\AbeComponentCheck\20200720\zernCoeffs_estimated_', ...
    num2str(zNum), '.txt']; % Zernike coefficients file
coeffsRaw1 = importdata(flieZernCoeffs1);
coeffsRaw2 = importdata(flieZernCoeffs2);
% zernRawNum = length(coeffsRaw);
% coeffsSigns = ones(1,zernRawNum); % to correct the mismatch between MATLAB and HASO
% coeffsSigns([2,5, 7, 10, 12, 14, 15, 17, 19, 21, 23, 26, 28, 30, 32]) = -1; % need further update
% coeffs = coeffsRaw .* coeffsSigns;

% coeffs = coeffs * 0;
% coeffs(8) = -1;
% Imagine Opitc convetion to Fringe convention
% 1st: tilt <--> 1st: piston
zTotal = 17;
pIn = 5:4+zTotal;
Sx = 512;
Sy = Sx;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
coeffs0 = zeros(1,17);
coeffs0(zNum-3) = 0.3;
coeffs1 = -coeffsRaw1(4:20)*0.6;
coeffs2 = -coeffsRaw2(4:20);
x = 4:3+zTotal;
figure,
plot(x,coeffs0,x,coeffs1,x,coeffs2,'LineWidth',2); % imagine optic convention
legend( 'ideal','measured','estimated');
xlabel('Zernike Coeff Index');
ylabel('Zernike Coeff Magnitude');
set(gca,'FontSize', 14);
title('Zernike coefficients');
% show wavefront
waveFront = nan(Sx, Sy);
waveFrontVector = create_wavefront(pIn,coeffs0,r(idx),theta(idx));
wMin = min(waveFrontVector(:));
wMax = max(waveFrontVector(:));
waveFront(idx) = waveFrontVector;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Ideal');
waveFront = nan(Sx, Sy);
waveFrontVector = create_wavefront(pIn,coeffs1,r(idx),theta(idx));
waveFront(idx) = waveFrontVector;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Measured');
waveFront = nan(Sx, Sy);
waveFrontVector = create_wavefront(pIn,coeffs2,r(idx),theta(idx));
waveFront(idx) = waveFrontVector;
xi = 1:Sx;
figure;
pcolor(xi,xi,waveFront), shading interp
axis square, Fc = colorbar;
caxis([wMin wMax])
xlabel(Fc,'\mum');
title('Estimated');

caxis([wMin wMax])
