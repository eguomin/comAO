% close all
zernNum = 6;
upLim = 0.1;
bottomLim = -0.1;
tNum = 300;
r0p5 = 0.5;
r1 = 1;
r2 = 2;
r3 = 3;
for i = 4:zernNum
    y0 = coeffsRaw_p5x1_1(1:tNum,i);
    y0 = y0(:) - mean(y0(:));
    y = y0;
    x0 = 1:tNum;
    x =  x0 * r0p5;
    figure,plot(x,y, 'LineWidth',2);
    y1 = reshape(y0,[2, tNum/2]);
    y = mean(y1,1);
    x0 = 1:length(y);
    x =  x0 * 1;
    hold on, plot(x,y, 'LineWidth',2);
    y1 = reshape(y0,[4, tNum/4]);
    y = mean(y1,1);
    x0 = 1:length(y);
    x =  x0 * 2;
    hold on, plot(x,y, 'LineWidth',2);
    y1 = reshape(y0,[6, tNum/6]);
    y = mean(y1,1);
    x0 = 1:length(y);
    x =  x0 * 3;
    hold on, plot(x,y, 'LineWidth',2);
    legend('0.5 s','1 s', '2 s', '3 s');
    xlabel('Time (s)');
    ylabel('Magnitude');
    xlim([0 150]);
    ylim([bottomLim upLim]);
end