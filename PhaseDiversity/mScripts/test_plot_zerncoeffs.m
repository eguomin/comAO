close all;
% import ground truth values: CheckList

for i = 1:32
% i = 1;
flieZernCoeffs = ['D:\Data\20200627_AbeComponentCheck\ZernIdx_', num2str(i),'_diff.txt']; % Zernike coefficients file
coeffsRaw = importdata(flieZernCoeffs);
zernRawNum = length(coeffsRaw);
x = 1:zernRawNum;
coeffsGT = zeros(size(coeffsRaw));
coeffsGT(i) = CheckList(i);
if(coeffsRaw(i)<0)
    coeffsRaw = -coeffsRaw;
end
% coeffsGT(i) = CheckListGT(i)*coeffsSigns(i);
% coeffsRaw = coeffsRaw/coeffsRaw(i);
figure,plot(x,coeffsRaw(:), 'LineWidth',2);
hold on, plot(x,coeffsGT(:), 'LineWidth',2);
upLim = max([max(coeffsRaw(:)) max(coeffsGT(:))]) * 1.1;
bottomLim = min([min(coeffsRaw(:)) min(coeffsGT(:))]) * 1.1;
xlim([0 33]);
ylim([bottomLim upLim]);
legend('measured', 'ground truth');
set(gca,'FontSize',18);
end