pIn = 1:20;
Sx = 256;
pixelSize = 0.096;
lambda = 0.532;
NA = 1.2;

fileZernIn = 'D:\Program\TempData\Comp_G_P\beads_31\results256\coef_gt_length.txt';
coeffsRaw = importdata(fileZernIn);
coeffsGT = zernidx2idxcoeffs(4:21, coeffsRaw, 'NOLL2ANSI', 20)';

fileZernIn = 'D:\Program\TempData\Comp_G_P\beads_31\results256\coef_est_g_length.txt';
coeffsRaw = importdata(fileZernIn);
coeffsG = zernidx2idxcoeffs(4:21, coeffsRaw, 'NOLL2ANSI', 20)';

fileZernIn = 'D:\Program\TempData\Comp_G_P\beads_31\results256\coef_est_p_length.txt';
coeffsRaw = importdata(fileZernIn);
coeffsP = zernidx2idxcoeffs(4:21, coeffsRaw, 'NOLL2ANSI', 20)';

fileZernIn = 'D:\Program\TempData\Comp_G_P\beads_31\beads_5ImgCorr_31\zernCoeffs_estimated.txt';
%coeffsRaw = importdata(fileZernIn);
coeffsRaw = fileIO_lvtxt2idx(fileZernIn);
coeffsM = -coeffsRaw(pIn)';

[~, ~, staGT, ~] = coeffs2wavefront(pIn,coeffsGT,Sx,pixelSize, lambda, NA, 0);
[~, ~, staG, ~] = coeffs2wavefront(pIn,coeffsGT - coeffsG,Sx,pixelSize, lambda, NA, 0);
[~, ~, staP, ~] = coeffs2wavefront(pIn,coeffsGT - coeffsP,Sx,pixelSize, lambda, NA, 0);
[~, ~, staM, ~] = coeffs2wavefront(pIn,coeffsGT - coeffsM,Sx,pixelSize, lambda, NA, 0);