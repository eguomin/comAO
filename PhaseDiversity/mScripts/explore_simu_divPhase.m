% explore the effect of type/amplitude of the diversity phase based on simulations

% By: Min Guo
% Mar. 10, 2020;
% Modifications: Aug. 19, 2020

clear all;
close all
% % % default settings
tStart = tic;
flagGPU = 0;
if(flagGPU)
    devNum = 1;
    gpuDevice(devNum);
end
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
pIn = 3:35; % 0: piston; 1:tilt X; 2: tilt Y;
nType = 'both';
SNR = 20;

% % % Unknown aberrations: to be estimated 
flagSetCoeffs = 1; % Set input Zernike coefficients:
if(flagSetCoeffs==1)
    coeffsFile = '..\..\..\computAO\Data\abeAndNoise\aberrationRondom2_0p15.mat';
    varNames = cell(1,3);
    varNames{1} = 'pIn';
    varNames{2} = 'aValue';
    varNames{3} = 'coeffs';
    importfile(coeffsFile,varNames);
    zernNum = length(pIn);
    expNum = size(coeffs,1);
elseif(flagSetCoeffs==2)
    expNum = 50;
    zernNum = length(pIn);
    coeffs = zeros(expNum,zernNum);
    aValue = 0.15;
    zernType = 'random2'; % 'defocus','astig', 'coma','trefoil','sphe', 'random1'
    for i = 1:expNum      
        % % % higher weight for lower order indices
        zernType = 'random2';
        c = gen_zern_coeffs(pIn,aValue,zernType);
        coeffs(i,:) = c;
    end
end

% % % diversity images:
imgNum = 2; % 2 or 3
abType1 = 'astig'; % 'astig', 'defocus', 'sphe', 'idx', or xx
abValue1 = sqrt(6); % defocus: Z3 = sqrt(3), astig: Z5 = sqrt(6)
            % coma: Z8 = sqrt(8), trefoil: Z9 = sqrt(8)
            % sphe: Z12 = sqrt(5)
            % idx: customize aberrations
flagSmoothEdge = 0;

% % % reconstruction settings
itLimit = 10; % max iteration number
gamma = 1e-10; % gamma range: 1e-14 ~ 1e-6
alpha = 0; % alpha: 0
penalChioce = 1;
% bootstrapping steps,e.g., 4th (14), 5th(20), 6th(27) and 7th(35) orders
zernSteps = [14 20 27 pIn(zernNum)] - pIn(1) + 1; 

% % % files and folders
fileFolderIn = '..\DataForTest\Simu\';
fileNameIn = 'imSample1';
fileImgSample = [fileFolderIn, fileNameIn, '.tif'];
img0 = single(ReadTifStack(fileImgSample));
[Sx, ~] = size(img0);
fileFolderOut = ['..\..\..\computAO\Data\exploreDivPhase\',...
    fileNameIn, '\'];
fileNameOut = [abType1,'_', nType, '_SNR_', num2str(SNR)];
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

% Known aberration: for diverisity image 1
if(strcmp(abType1,'idx'))
    zernIdx = [4 5];
    zernValue = [0.2 0.5];
    coeffs_temp = gen_zern_coeffs(pIn,abValue1,abType1,zernIdx,zernValue);
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,coeffs_temp,Sx,...
            pixelSize, lambda, NA, 0);
    coeffs_delta = coeffs_temp/staPara.rmsLength;
else
    coeffs_delta = gen_zern_coeffs(pIn,abValue1,abType1);
end
pdNum = size(coeffs_delta, 1);

% deversity phase amount settings
deltaRMSs = [0.05 0.1:0.1:1.1 1.2:0.2:2];
rNum = length(deltaRMSs);

% simulation
pIn0 = pIn;
coeffsIns = zeros(expNum,zernNum, rNum);
coeffsEsts = zeros(expNum,zernNum, rNum);
coeffs_deltas = zeros(pdNum,zernNum, rNum);
staWaveFrontIn = zeros(expNum,2, rNum);
staWaveFrontEst = zeros(expNum,2, rNum);
staWaveFrontEstError = zeros(expNum,2, rNum);
for j = 1:rNum 
    % cTime1 = toc(tStart);
    disp('********************************************************');
    disp(['Batch simulation #: ', num2str(j)]);
    coeffsIn0 = coeffs;
    coeffs_delta0 = coeffs_delta*(deltaRMSs(j));
    coeffsIns(:,:,j) = coeffsIn0;
    coeffs_deltas(:,:,j) = coeffs_delta0;
    gamma0 = gamma;
    alpha0 = alpha;
    coeffsEst = recon_zern_simuBatch(pIn0,coeffsIn0, img0, coeffs_delta0, ...
        gamma0, alpha0, itLimit, zernSteps, pixelSize, lambda, NA, nType, ...
        SNR, flagGPU, penalChioce);
    coeffsEsts(:,:,j) = coeffsEst;
    
    % % % statistics
    disp('Performing statistics...');
    for i = 1:expNum
        c = coeffs(i,:);
        cEst = coeffsEst(i,:);
        [~, ~, staPara, ~] = coeffs2wavefront(pIn,c,Sx,...
            pixelSize, lambda, NA, 0);
        staWaveFrontIn(i,1,j) = staPara.rmsLength;
        staWaveFrontIn(i,2,j) = staPara.pvLength;
        [~, ~, staPara, ~] = coeffs2wavefront(pIn,cEst,Sx,...
            pixelSize, lambda, NA, 0);
        staWaveFrontEst(i,1,j) = staPara.rmsLength;
        staWaveFrontEst(i,2,j) = staPara.pvLength;
        [~, ~, staPara, ~] = coeffs2wavefront(pIn,cEst-c,Sx,...
            pixelSize, lambda, NA, 0);
        staWaveFrontEstError(i,1,j) = staPara.rmsLength;
        staWaveFrontEstError(i,2,j) = staPara.pvLength;
    end
end
staWaveFrontIn_Ave = mean(staWaveFrontIn, 1);
staWaveFrontIn_SD = sqrt(var(staWaveFrontIn, 1));
staWaveFrontEst_Ave = mean(staWaveFrontEst, 1);
staWaveFrontEst_SD = sqrt(var(staWaveFrontEst, 1));
staWaveFrontEstError_Ave = mean(staWaveFrontEstError, 1);
staWaveFrontEstError_SD = sqrt(var(staWaveFrontEstError, 1));
cTime = toc(tStart);
save([fileFolderOut, fileNameOut,'.mat']);
disp(['Processing completed!!! Total time cost:', num2str(cTime), ' s']);
% show RMS and PV curves 
x = deltaRMSs;
y = staWaveFrontEstError_Ave(1,1,:);
figure, plot(x, y(:), 'LineWidth',2);
xlabel('RMS of Diversity Phase');
ylabel('RMS Error ({\mu}m)');
legend([nType, ' SNR=', num2str(SNR)]);
set(gca,'fontsize', 18);
savefig([fileFolderOut, fileNameOut, '_RMSCurve.fig']);
y = staWaveFrontEstError_Ave(1,2,:);
figure, plot(x, y(:), 'LineWidth',2);
xlabel('RMS of Diversity Phase');
ylabel('PV Error ({\mu}m)');
legend([nType, ' SNR=', num2str(SNR)]);
set(gca,'fontsize', 18);
savefig([fileFolderOut, fileNameOut, '_PVCurve.fig']);
