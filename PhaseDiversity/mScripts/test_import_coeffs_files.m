% test import coeffs files
pathIn = 'D:\Data\20200718_monitorBead\WF_processed\WF_abeFree_p5x1_5_Coeffs\';

fileNum = 1000;
coeffsRaw = zeros(fileNum, 32);
for i = 0:fileNum-1
    flieZernCoeffs = [pathIn, 'WF_', num2str(i),'.txt']; % Zernike coefficients file
    coeffs = importdata(flieZernCoeffs);
    coeffsRaw(i+1,:) = coeffs;
end

% rename coeffsRaw