% show zernike coefficients from imagine optics
% 
% 
close all;
zernNum = 20;
x = 1:zernNum;
pathMain = 'C:\Programs\computAO\Data\20200816_DMresponseAnalysis\';
r = 0.3;
folder1 = 'Batch_0627AOMI_1um_1um_0p3';
folder2 = 'Batch_0807AOMI_0p5um_1um_0p3';
folder3 = 'Batch_0816AOMI_1um_1um_0p3';
folder4 = 'Batch_0816_1sx2_0p5um_0p5um_0p3';
subFolder = 'ANSI_idx_';
fileName = 'Net_0p3_zernCoefNet.txt';
for i = 1:20
    zern0 = zeros(1,zernNum);
    zern0(i) = r;
    y = zern0;
    figure,plot(x,y, '-k', 'LineWidth',2);
    flieZernCoeffs = [pathMain, folder1, '\', subFolder, num2str(i), '\', fileName]; 
    coeffsRaw = importdata(flieZernCoeffs);
    y = -coeffsRaw(1:zernNum);
    hold on, plot(x,y, 'LineWidth',2);
    flieZernCoeffs = [pathMain, folder2, '\', subFolder, num2str(i), '\', fileName]; 
    coeffsRaw = importdata(flieZernCoeffs);
    y = -coeffsRaw(1:zernNum);
%     hold on, plot(x,y, 'LineWidth',1);
%     flieZernCoeffs = [pathMain, folder3, '\', subFolder, num2str(i), '\', fileName]; 
%     coeffsRaw = importdata(flieZernCoeffs);
%     y = -coeffsRaw(1:zernNum);
%     hold on, plot(x,y, 'LineWidth',1);
    flieZernCoeffs = [pathMain, folder4, '\', subFolder, num2str(i), '\', fileName]; 
    coeffsRaw = importdata(flieZernCoeffs);
    y = -coeffsRaw(1:zernNum);
    hold on, plot(x,y, 'LineWidth',2);
%     legend('input', '0627 1um 1um', '0807 0p5um 1um', '0816 1um 1um', '0816 0p5um 0p5um');
%     legend('input', '1um bead', '0.5um bead');
    xlabel('Zernike Index (ANSI)');
    ylabel('Coef Value');
    xlim([0 21]);
    ylim([-0.2 r+0.1]);
    set(gca,'fontsize', 18);
end