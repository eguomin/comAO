x = staWaveFrontIn_Ave(1,1,:);
y = staWaveFrontEstError_Ave(1,1,:);
% figure, plot(x(:), y(:), 'LineWidth',2);
hold on, plot(x(:), y(:), 'LineWidth',2);
xlabel('RMS of Aberration');
ylabel('RMS Error ({\mu}m)');
% legend([nType, ' SNR=', num2str(SNR), ' RMS=', num2str(abValue1)]);

legend('RMS=0.5','RMS=0.3','RMS=0.3+0.3', 'RMS=0.3, SNR=20');