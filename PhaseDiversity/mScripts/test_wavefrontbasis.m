% generate wavefront basis: ANSI 1-35
pIn = 1:35;
pixelSize = 0.096; % um
lambda = 0.550; % um
NA = 1.2;
Sx = 512;
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
waveFront = zeros(Sx,Sx);
waveFronts = zeros(Sx,Sx,length(pIn(:)));
conType = 'ANSI';
nFlag = 'none';
for i = 1:length(pIn(:))
    p = pIn(i);
    waveFront(idx) = create_wavefront(p, 1, r(idx), theta(idx), nFlag, conType);
    waveFronts(:,:,i) = waveFront;
end
WriteTifStack(waveFronts,'waveFrontBasis_ANSI.tif',32);