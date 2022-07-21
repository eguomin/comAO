% explore image quality metric IQM vs. RMS

% June 19, 2020
% clear all;
% close all
% parameters
tic;
pIn = 4:21; % 1: piston; 2:tilt X; 3: tilt Y;
zernNum = length(pIn);
zernCoeffOrder = max(pIn);
% nType = 'both';
nType = 'none';
SNR = 50;
expNum = 10;
flagNewAberration = 1;
amRatios = [0:0.2:0.8 1:1:5 8:3:25];
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
imgName = 'imSample1';
fileFolderOut = ['C:\Programs\computAO\Data\IQM\IQM_results\', imgName, '_', nType, num2str(SNR),'\'];
fileImgSample = ['C:\Programs\computAO\Data\IQM\', imgName, '.tif'];
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
% mType = 7;
% IQM_DCTSs = zeros(rNum, expNum);
mTypeNum = 7;
IQM_all = zeros(mTypeNum, rNum, expNum);
RMS_Values = zeros(rNum, expNum);

for j = 1:rNum 
    amR = amRatios(j);
    coeffsIn1 = coeffsIn * amR;
    for i = 1: expNum
        coeffsIn0 = coeffsIn1(i,:);
        [img, ~, ~] = gen_simu_images(img0, pIn, coeffsIn0, pixelSize, lambda, NA, nType, SNR);
        if(~strcmp(nType, 'none'))
            img = addnoiseSNR(img,nType,SNR);
        end
        for k = 1:mTypeNum
            IQM_all(k,j,i) = calculate_iqm(img,k, 3);
        end
        [~, ~, staPara, ~] = coeffs2wavefront(pIn,coeffsIn0,Sx, pixelSize, lambda, NA, 0);
        RMS_Values(j,i) = staPara.rmsLength;
    end   
end
ave_IQM_all = mean(IQM_all, 3);
ave_RMS_Values = mean(RMS_Values, 2);
cTime = toc
% plot 1, 3, 4, 6, 7
figure;
for nPlot = [1, 3, 4, 6, 7]
    x = ave_RMS_Values(:);
    y = ave_IQM_all(nPlot,:);
    y = y(:) - min(y(:));
    y = y(:)/max(y(:));
    hold on, plot(x, y, 'LineWidth',4);
end
% legend('Max Intensity', 'Contrast', 'Sharpness Image', 'Sharpness Fre', 'DCTS');
ylim([0 1.2]);
xlim([0 1.1]);
set(gca,'FontSize',18)
save([fileFolderOut 'data.mat']);