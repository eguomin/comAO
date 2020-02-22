% explore amplitude of the aberratin/diversity phase based on simulations
% clear all;
% close all
pathOutMain = 'I:\AO_data\20200220_Simu\divesityApm\hybrid_asti_coma\';
flagNewAberration = 0;
% parameters
pIn = 4:21; % 1: piston; 2:tilt X; 3: tilt Y;
zernNum = length(pIn);
iteration = 5; % note: more zernike orders --> more iterations? *******
zernCoeffOrder = max(pIn);
% nType = 'none';
% gamma = 1e-14;
nType = 'poisson';
gamma = 1e-6;
fileFolderOut = [pathOutMain , nType, '\'];

expNum = 50;
pdNum = 1;
% amRatios = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8 10];
amRatios = [0.2:0.2: 0.8 1:1:10];
rNum = length(amRatios);
% random zernike coefficients
cAm = 0.3;
if(flagNewAberration)
    coeffsIn = zeros(expNum,zernNum);
    for i = 1:expNum
        cRandom = cAm.*rand(1,zernNum) - cAm/2;
        coeffsIn(i,:) = cRandom;
    end
end
% diversity phase
dAm = 0.3;
coeffs_deltaInitial = zeros(4,zernNum);
% coeffs_deltaInitial(1,5-pIn(1)+1) = dAm; % astigmitism;
coeffs_deltaInitial(1,5-pIn(1)+1) = dAm; % hybrid;
coeffs_deltaInitial(1,7-pIn(1)+1) = dAm; % hybrid;
coeffs_deltaInitial(2,7-pIn(1)+1) = dAm; % coma;
coeffs_deltaInitial(3,4-pIn(1)+1) = dAm; % defocus;
coeffs_deltaInitial(4,9-pIn(1)+1) = dAm; % spherical;
% other parameters
pixelSize = 0.096; % um
lambda = 0.550; % um
NA = 1.2;
flagGPU = 1;
fileImgSample = '..\DataForTest\Simu\imSample2.tif';
img0 = double(ReadTifStack(fileImgSample));
[Sx,Sy] = size(img0);
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

% simulation
pIn0 = pIn;
coeffs_delta = coeffs_deltaInitial(1:pdNum,:);
coeffsIns = zeros(expNum,zernNum, rNum);
coeffsOuts = zeros(expNum,zernNum, rNum);
coeffs_deltas = zeros(pdNum,zernNum, rNum);
% coeffsIn = coeffsIn;
for j = 1:rNum 
    amR = amRatios(j);
    coeffsIn0 = coeffsIn;
    coeffs_delta0 = coeffs_delta.*amR;
    coeffsIns(:,:,j) = coeffsIn0;
    coeffs_deltas(:,:,j) = coeffs_delta0;
    [coeffsOut, cTime] = recon_zern_simuBatch(pIn0,coeffsIn0, img0, coeffs_delta0, ...
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
rmsDelta = zeros(rNum,1);
pvDelta = zeros(rNum,1);
coeffsErrors = coeffsOuts - coeffsIns;
for j = 1:rNum
    coeffsIn0 = coeffsIns(:,:,j);
    coeffsOut0 = coeffsOuts(:,:,j);
    coeffsError0 = coeffsOut0 - coeffsIn0;
    coeffs_delta0 = coeffs_deltas(:,:,j);
    for i = x
        [~, ~, staParaIn, ~] = coeffs2wavefront(pIn0,coeffsIn0(i,:),Sx, pixelSize, lambda, NA, 0);
        [~, ~, staParaOut, ~] = coeffs2wavefront(pIn0,coeffsOut0(i,:),Sx, pixelSize, lambda, NA, 0);
        [~, ~, staParaError, ~] = coeffs2wavefront(pIn0,coeffsError0(i,:),Sx, pixelSize, lambda, NA, 0);
        rmsIn(j,i) = staParaIn.rmsLength;
        pvIn(j,i) = staParaIn.pvLength;
        rmsOut(j,i) = staParaOut.rmsLength;
        pvOut(j,i) = staParaOut.pvLength;
        rmsError(j,i) = staParaError.rmsLength;
        pvError(j,i) = staParaError.pvLength;
    end
    [~, ~, staParaDelat, ~] = coeffs2wavefront(pIn0,coeffs_delta0(1,:),Sx, pixelSize, lambda, NA, 0);
    rmsDelta(j,1) = staParaDelat.rmsLength;
    pvDelta(j,1) = staParaDelat.pvLength;
end
rmsInAve = mean(rmsIn, 2);
pvInAve = mean(pvIn, 2);
rmsOutAve = mean(rmsOut,2);
pvOutAve = mean(pvOut,2);
rmsErrorAve = mean(rmsError,2);
pvErrorAve = mean(pvError,2);
save([fileFolderOut 'data.mat']);

figure, plot(rmsInAve(:),rmsErrorAve(:),rmsDelta(:),rmsErrorAve(:));
title('Measure in RMS (\mum)');
legend('Aberrated','Diveristy');
xlabel('RMS of diversity or aberrated');
ylabel('RMS Error');
savefig([fileFolderOut 'RMScurve.fig']);
figure, plot(pvInAve(:),pvErrorAve(:),pvDelta(:),pvErrorAve(:));
title('Measure in PV (\mum)');
legend('Aberrated','Diveristy');
xlabel('PV of diversity or aberrated');
ylabel('PV Error');
savefig([fileFolderOut 'PVcurve.fig']);
% show typical results at the best reconstruction
[y, idx] = min(rmsErrorAve(:));
% coeffsInAve = mean(coeffsIns(:,:,idx), 1);
% coeffsOutAve = mean(coeffsOuts(:,:,idx), 1);
% coeffsErrorAve = coeffsOutAve - coeffsInAve;
% coeffsDeltaAve = mean(coeffs_deltas(:,:,j), 1);
coeffsInAve = coeffsIns(1,:,idx);
coeffsOutAve = coeffsOuts(1,:,idx);
coeffsErrorAve = coeffsOutAve - coeffsInAve;
coeffsDeltaAve = coeffs_deltas(1,:,j);
[waveFrontIn, ~, ~, ~] = coeffs2wavefront(pIn0,coeffsInAve(:),Sx, pixelSize, lambda, NA, 0);
[waveFrontOut, ~, ~, ~] = coeffs2wavefront(pIn0,coeffsOutAve(:),Sx, pixelSize, lambda, NA, 0);
[waveFrontError, ~, ~, ~] = coeffs2wavefront(pIn0,coeffsErrorAve(:),Sx, pixelSize, lambda, NA, 0);
[waveFrontDelta, ~, ~, ~] = coeffs2wavefront(pIn0,coeffsDeltaAve(:),Sx, pixelSize, lambda, NA, 0);
xi = 1: Sx;
figure, subplot(2,2,1);
pcolor(xi,xi,waveFrontIn), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: aberrated');
subplot(2,2,2);
pcolor(xi,xi,waveFrontDelta), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: diversity');
subplot(2,2,3);
pcolor(xi,xi,waveFrontOut), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: estimated');
subplot(2,2,4);
pcolor(xi,xi,waveFrontError), shading interp
axis square, Fc = colorbar;
xlabel(Fc,'\mum');
title('Wavefront: error');
savefig([fileFolderOut 'wavefrontExample.fig']);
% % % % % % 
% figure, plot(x,rmsIn,x,rmsOut(1,:)); 
% title('RMS (\mum)');
% for j = 2:rNum
%     hold on, plot(x,rmsOut(j,:));
% end
% legend('0','0.5','1','2','5', '10');
% 
% figure, plot(x,pvIn,x,pvOut(1,:)); 
% title('PV (\mum)');
% for j = 2:rNum
%     hold on, plot(x,pvOut(j,:));
% end
% legend('0','0.5','1','2','5', '10');
% figure, plot(x,rmsError(1,:))
% for j = 2:rNum
%     hold on, plot(x,rmsError(j,:));
% end
% title('RMS Error (\mum)');
% legend('0.5','1','2','5', '10');
% figure,plot(x,pvError(1,:)); 
% for j = 2:rNum
%     hold on, plot(x,pvError(j,:));
% end
% title('PV Error (\mum)');
% legend('0.5','1','2','5', '10');