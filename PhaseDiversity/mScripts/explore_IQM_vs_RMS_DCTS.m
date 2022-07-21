% explore image quality metric DCTS vs. RMS

% June 4, 2020
% clear all;
% close all
% parameters
pIn = 4:21; % 1: piston; 2:tilt X; 3: tilt Y;
zernNum = length(pIn);
zernCoeffOrder = max(pIn);
nType = 'poisson';
SNR = 2;
expNum = 10;
flagNewAberration = 1;
amRatios = [0:0.2:0.8 1:1:10 12:2:20];
rNum = length(amRatios);

% generate zernike coefficients
cAm = 0.1;
if(flagNewAberration)
    coeffsIn = zeros(expNum,zernNum);
    for i = 1:expNum
        cRandom = cAm.*rand(1,zernNum) - cAm/2;
        coeffsIn(i,:) = cRandom;
    end
end

% input images and output folder
fileFolderOut = 'C:\Programs\computAO\Data\IQM\DCTS\imSample1_SNR2\';
fileImgSample = 'C:\Programs\computAO\Data\IQM\imSample1.tif';
img0 = double(ReadTifStack(fileImgSample));
[Sx,Sy] = size(img0);
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

% Simulate images: calculate DCTS and RMS 
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
mType = 7;
IQM_DCTSs = zeros(rNum, expNum);
RMS_Values = zeros(rNum, expNum);
for j = 1:rNum 
    amR = amRatios(j);
    coeffsIn1 = coeffsIn * amR;
    for i = 1: expNum
        coeffsIn0 = coeffsIn1(i,:);
        [img, ~, ~] = gen_simu_images(img0, pIn, coeffsIn0, pixelSize, lambda, NA, nType);
        if(~strcmp(nType, 'none'))
            img = addnoiseSNR(img,nType,SNR);
        end
        IQM_DCTSs(j,i) = calculate_iqm(img,mType);
        [~, ~, staPara, ~] = coeffs2wavefront(pIn,coeffsIn0,Sx, pixelSize, lambda, NA, 0);
        RMS_Values(j,i) = staPara.rmsLength;
    end   
end
ave_IQM_DCTS = mean(IQM_DCTSs, 2);
ave_RMS_Values = mean(RMS_Values, 2);
figure, plot(ave_RMS_Values(:), ave_IQM_DCTS(:));
save([fileFolderOut 'data.mat']);