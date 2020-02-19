% explore amplitude of the diversity phase based on simulations
clear all
% parameters
pIn = 4:21; % 1: piston; 2:tilt X; 3: tilt Y;
zernNum = length(pIn);
iteration = 5; % note: more zernike orders --> more iterations? *******
zernCoeffOrder = max(pIn);
nType = 'poisson';
gamma = 1e-8;

expNum = 20;
pdNum = 2;
% amRatios = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8 10];
amRatios = [0.5, 1, 2, 5, 10];
rNum = length(amRatios);
% random zernike coefficients
cAm = 0.3;
coeffsIn = zeros(expNum,zernNum);
for i = 1:expNum
    cRandom = cAm.*rand(1,zernNum) - cAm/2;
    coeffsIn(i,:) = cRandom;
end
% diversity phase
dAm = 0.3;
coeffs_delta = zeros(4,zernNum);
coeffs_delta(1,5-pIn(1)+1) = dAm; % astigmitism;
coeffs_delta(2,7-pIn(1)+1) = dAm; % coma;
coeffs_delta(3,4-pIn(1)+1) = dAm; % defocus;
coeffs_delta(4,9-pIn(1)+1) = dAm; % spherical;
% other parameters
pixelSize = 0.096; % um
lambda = 0.550; % um
NA = 1.2;
flagGPU = 1;
fileImgSample = '..\DataForTest\Simu\imSample2.tif';
img0 = double(ReadTifStack(fileImgSample));
[Sx,Sy] = size(img0);

% simulation
pIn0 = pIn;
coeffs_delta0 = coeffs_delta(1:pdNum,:);
coeffsOuts = zeros(expNum,zernNum, rNum);
for j = 1:rNum 
    amR = amRatios(j);
    coeffsIn0 = coeffsIn;
    [coeffsOut, cTime] = recon_zern_simuBatch(pIn0,coeffsIn0, img0, amR.*coeffs_delta0, ...
        zernCoeffOrder, iteration, gamma, pixelSize, lambda, NA, nType, flagGPU);
    coeffsOuts(:,:,j) = coeffsOut;
end
% analysis
% figure, plot(pIn, coeffsIn(1,:),pIn, coeffsOut(1,:));
% legend('GroundTruth', 'Estimated');
% title('Zernike coeffs (first case) (\mum)');
[expNum, ~] = size(coeffsIn0);
x = 1:expNum;
rmsIn = zeros(1,expNum);
pvIn = zeros(1,expNum);
rmsOut = zeros(rNum,expNum);
pvOut = zeros(rNum,expNum);
rmsError = zeros(rNum,expNum);
pvError = zeros(rNum,expNum);
for j = 1:rNum
    coeffsOut = coeffsOuts(:,:,j);
    coeffsError = coeffsOut - coeffsIn0;
    for i = x
        [~, ~, staParaIn, ~] = coeffs2wavefront(pIn0,coeffsIn0(i,:),Sx, pixelSize, lambda, NA, 0);
        [~, ~, staParaOut, ~] = coeffs2wavefront(pIn0,coeffsOut(i,:),Sx, pixelSize, lambda, NA, 0);
        [~, ~, staParaError, ~] = coeffs2wavefront(pIn0,coeffsError(i,:),Sx, pixelSize, lambda, NA, 0);
        rmsIn(1,i) = staParaIn.rmsLength;
        pvIn(1,i) = staParaIn.pvLength;
        rmsOut(j,i) = staParaOut.rmsLength;
        pvOut(j,i) = staParaOut.pvLength;
        rmsError(j,i) = staParaError.rmsLength;
        pvError(j,i) = staParaError.pvLength;
    end
end

figure, plot(x,rmsIn,x,rmsOut(1,:)); 
title('RMS (\mum)');
for j = 2:rNum
    hold on, plot(x,rmsOut(j,:));
end
legend('0','0.5','1','2','5', '10');

figure, plot(x,pvIn,x,pvOut(1,:)); 
title('PV (\mum)');
for j = 2:rNum
    hold on, plot(x,pvOut(j,:));
end
legend('0','0.5','1','2','5', '10');
figure, plot(x,rmsError(1,:))
for j = 2:rNum
    hold on, plot(x,rmsError(j,:));
end
title('RMS Error (\mum)');
legend('0.5','1','2','5', '10');
figure,plot(x,pvError(1,:)); 
for j = 2:rNum
    hold on, plot(x,pvError(j,:));
end
title('PV Error (\mum)');
legend('0.5','1','2','5', '10');