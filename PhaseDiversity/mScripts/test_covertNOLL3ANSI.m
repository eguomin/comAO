fileZernIn = 'D:\Program\TempData\Comp_G_P\beads_2\results\coef_gt_length.txt';
coeffsRaw = importdata(fileZernIn);
coeffsANSI = zernidx2idxcoeffs(4:21, coeffsRaw, 'NOLL2ANSI', 20)';

