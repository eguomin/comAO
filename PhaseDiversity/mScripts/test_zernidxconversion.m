% test zernike index conversion
p0 = 1:50;
[n,m] = zernidx2nm(p0,'ANSI');
p = zernidx2idx(p0,'ANSI2Wyant');
p0nmpANSI2Wyant = [p0;n;m;p;]';
[n,m] = zernidx2nm(p0,'Wyant');
p = zernidx2idx(p0,'Wyant2ANSI');
p0nmpWyant2ANSI = [p0;n;m;p;]';

coeffsANSI2Wyant2 = zernidx2idxcoeffs(1:28, 1:28, 'ANSI2Wyant', 32)';
coeffsWyant2ANSI2 = zernidx2idxcoeffs(1:32, 1:32, 'Wyant2ANSI', 28)';