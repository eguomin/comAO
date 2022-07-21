zernNum = 32;
upLim = 0.1;
bottomLim = -0.1;
x1 = 1:300;
x2 = 1:150;
r0p5 = 0.5;
r1 = 1;
r2 = 2;
r3 = 3;
for i = 1:zernNum
    x = x1 * r0p5;
    y = coeffsRaw_p5(:,i);
    y = y(:) - mean(y(:));
    figure,plot(x,y, 'LineWidth',2);
    x = x1 * r1;
    y = coeffsRaw_1(:,i);
    y = y(:) - mean(y(1:150));
    hold on, plot(x,y, 'LineWidth',2);
    x = x2 * r2;
    y = coeffsRaw_2(:,i);
    y = y(:) - mean(y(1:75));
    hold on, plot(x,y, 'LineWidth',2);
    x = x2 * r3;
    y = coeffsRaw_3(:,i);
    y = y(:) - mean(y(1:50));
    hold on, plot(x,y, 'LineWidth',2);
    legend('0.5 s','1 s', '2 s', '3 s');
    xlim([0 150]);
    ylim([bottomLim upLim]);
end

for i = 1:zernNum
    x = x1 * r0p5;
    y = coeffsRaw_p5_2(:,i);
    y = y(:) - mean(y(:));
    figure,plot(x,y, 'LineWidth',2);
    x = x1 * r1;
    y = coeffsRaw_1_2(:,i);
    y = y(:) - mean(y(1:150));
    hold on, plot(x,y, 'LineWidth',2);
    x = x2 * r2;
    y = coeffsRaw_2_2(:,i);
    y = y(:) - mean(y(1:75));
    hold on, plot(x,y, 'LineWidth',2);
    x = x2 * r3;
    y = coeffsRaw_3_2(:,i);
    y = y(:) - mean(y(1:50));
    hold on, plot(x,y, 'LineWidth',2);
    legend('0.5 s','1 s', '2 s', '3 s');
    xlim([0 150]);
    ylim([bottomLim upLim]);
end