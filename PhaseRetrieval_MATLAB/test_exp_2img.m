% test recon_zern_fileIO
clear all;
% close all;
fileFolderMain = 'D:\Data\20200212_PR\';
fileName = 'beads_3';
fileFolderIn = fileFolderMain;
imgNum = 2;
fileFolderOut = ['D:\Data\20200212_PR\', fileName, '_', num2str(imgNum), '_new\'];
repNum = 2;
zernCoeffOrder = 15;
iteration = 10;
gamma = 1e-6;
cropSize = 384;
bgValue = 300;
flagShowInput = 0;
flagShowRecon = 1;
flagGPU = 1;
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
[coeffsOut, cRMSE, cTime] = recon_zern_fileIO(fileFolderOut,fileFolderIn, fileName, imgNum, repNum,...
    zernCoeffOrder, iteration, gamma, cropSize, bgValue,flagShowInput, flagShowRecon, flagGPU);
