% A simple example: test phase diveristy based on experimental datasets.
%   1) recon_zern_fileIO function takes care of file IO and wavefront estimation.
%   2) estimated wavefronts (Zernike coefficients) are sent to deformable mirror 
%      to compensate the aberration.
%
% Current setup has been verified for the exp datasets 'bead_2' and 'bead_8'

% By: Min Guo, ShroffLab@NIH
% Last update: Mar. 11, 2022

clear all;
close all;
fileFolderMain = '..\DataForTest\Exp\';
fileName = 'beads_2';
fileFolderIn = [fileFolderMain, fileName, '\'];
imgNum = 2;
fileFolderOut = [fileFolderIn, 'TestResult_', num2str(imgNum), 'img\'];
repNum = 1;
zernIdxLimit = 20;
itLimit = 10;
gamma = 1e-10;
alpha = 0;
cropSize = 512;
bgValue = 300;
flagShowInput = 0;
flagShowRecon = 1;
flagGPU = 0;
idxType ='ANSI';
penalChioce = 1;
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
[coeffsOut, cRMSE, cTime] = recon_zern_fileIO(fileFolderOut,fileFolderIn, fileName, imgNum, repNum,...
    zernIdxLimit, itLimit, gamma, alpha, cropSize, bgValue,flagShowInput, flagShowRecon, flagGPU,...
    idxType,penalChioce);
