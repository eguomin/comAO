% test recon_zern_fileIO
clear all;
% close all;
fileFolderMain = '..\DataForTest\Exp\';
fileName = 'beads_25';
fileFolderIn = [fileFolderMain, fileName, '\'];
imgNum = 2;
fileFolderOut = [fileFolderIn, 'TestResult_', num2str(imgNum), 'img\'];
repNum = 2;
zernCoeffOrder = 15;
iteration = 10;
gamma = 1e-6;
cropSize = 384;
bgValue = 290;
flagShowInput = 0;
flagShowRecon = 1;
flagGPU = 0;
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
[coeffsOut, cRMSE, cTime] = recon_zern_fileIO(fileFolderOut,fileFolderIn, fileName, imgNum, repNum,...
    zernCoeffOrder, iteration, gamma, cropSize, bgValue,flagShowInput, flagShowRecon, flagGPU);
