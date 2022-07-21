% test_addnoise
% June 19, 2020

% nType = 'gaussian';
% nType = 'poisson';
nType = 'both';
% nType = 'none';
SNR = 20;
G = 10;
sType = 3;
% input images and output folder
imgName = 'imSample2';
fileFolderOut = ['C:\Programs\computAO\Data\testNoise\', imgName,'\'];
fileImgSample = ['C:\Programs\computAO\Data\groundTruth\', imgName, '.tif'];
img0 = double(ReadTifStack(fileImgSample));
[Sx,Sy] = size(img0);
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
img = addnoiseSNR(img0, nType, SNR, G, sType);
WriteTifStack(img, [fileFolderOut, nType, '_', num2str(SNR), '.tif'], 32);