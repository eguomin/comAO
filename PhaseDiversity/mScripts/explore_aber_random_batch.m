% explore generating random aberrations

expNum = 50;
pIn = 3:35;
zernNum = length(pIn);
coeffs = zeros(expNum,zernNum);
staWaveFront = zeros(expNum,2); 
Sx = 512;
pixelSize = 0.096; % um
lambda = 0.532; % um
NA = 1.2;
fileName = 'aberrationRondom2_0p10';
fileFolderOut = '..\..\..\computAO\Data\abeAndNoise\';
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

aValue = 0.10;
for i = 1:expNum
    % % % single component
    % zernType = 'defocus';
    % c = gen_zern_coeffs(pIn,aValue,zernType);
    
    % % % uniform weight for all indices
    % zernType = 'random1';
    % c = gen_zern_coeffs(pIn,aValue,zernType);
    
    % % % higher weight for lower order indices
    zernType = 'random2';
    c = gen_zern_coeffs(pIn,aValue,zernType);

    % % % statistics
    [~, ~, staPara, ~] = coeffs2wavefront(pIn,c,Sx,...
    pixelSize, lambda, NA, 0); 
    coeffs(i,:) = c;
    staWaveFront(i,1) = staPara.rmsLength;
    staWaveFront(i,2) = staPara.pvLength;
end
staWaveFront_Ave = mean(staWaveFront, 1);
staWaveFront_SD = sqrt(var(staWaveFront, 1));
save([fileFolderOut, fileName, '.mat']);