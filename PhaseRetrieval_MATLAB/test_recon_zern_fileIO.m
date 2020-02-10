% test recon_zern_fileIO

fileFolderOut = 'C:\Users\eguom\Documents\Test\';
fileFolderIn = 'C:\Users\eguom\Documents\GitHub\computationalAO\DataForTest\Exp\';
fileName = 'beads_1';
imgNum = 2;
repNum = 2;
zernCoeffOrder = 15;
iteration = 10;
gamma = 1e-6;
cropSize = 384;
bgValue = 290;
flagShowInput = 0;
flagShowRecon = 1;
flagGPU = 1;
[coeffsOut, cTime] = recon_zern_fileIO(fileFolderOut,fileFolderIn, fileName, imgNum, repNum,...
    zernCoeffOrder, iteration, gamma, cropSize, bgValue,flagShowInput, flagShowRecon, flagGPU);
cTime
