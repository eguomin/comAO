% test adding aberrations and noise to images

pIn = 3:35;
zernNum = length(pIn);
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
% % % noise settings
nType = 'both'; % 'gaussian', 'poisson'
% nType = 'none';
SNR = 20;
G = 10;
sType = 3;
% % % aberration
aValue = 0.2;
zernType = 'defocus'; % 'defocus','astig', 'coma','trefoil','sphe', 'random1'
coeffs = gen_zern_coeffs(pIn,aValue,zernType);
% % % input images and output folder
imgName = 'imSample1';
fileFolderOut = ['..\..\..\computAO\Data\abeAndNoise\', imgName,'\'];
fileImgSample = ['..\DataForTest\Simu\', imgName, '.tif'];
img0 = double(ReadTifStack(fileImgSample));
[Sx,Sy] = size(img0);
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

[~, ~, staPara, ~] = coeffs2wavefront(pIn,coeffs,Sx,...
    pixelSize, lambda, NA, 0);
WF_rms = staPara.rmsLength;
WF_pv = staPara.pvLength;
WriteTifStack(img0, [fileFolderOut, 'groundtruth.tif'], 32);
[imgs, ~, ~] = gen_simu_images(img0, pIn, zeros(1,zernNum), pixelSize, lambda, NA, 'none',SNR);
WriteTifStack(imgs, [fileFolderOut, 'none_abeFree.tif'], 32);
[imgs, ~, ~] = gen_simu_images(img0, pIn, coeffs, pixelSize, lambda, NA, 'none',SNR);
WriteTifStack(imgs, [fileFolderOut, 'none_', zernType, '_', num2str(aValue), '.tif'], 32);
[imgs, ~, ~] = gen_simu_images(img0, pIn, coeffs, pixelSize, lambda, NA, nType,SNR);
WriteTifStack(imgs, [fileFolderOut, nType, '_', num2str(SNR),... 
    '_', zernType, '_', num2str(aValue), '.tif'], 32);
save([fileFolderOut, nType, '_', num2str(SNR),... 
    '_', zernType, '_', num2str(aValue), '.mat']);