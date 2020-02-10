% create sample

% length unit
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

% Image and pixel size
Sx = 512;
Sy = 512;
pixelSize = 0.1625*um;
imSample = zeros(Sx, Sy);

subSx = 100;
subSy = subSx;

% point
% subImSample = zeros(subSx, subSy);
% subImSample(round((subSx+1)/2), round((subSy+1)/2)) = 1;
% subImSample1 = subImSample;

% points
iDistance = 20;
subImSample = zeros(subSx, subSy);
subImSample(round((subSx+1)/2)-iDistance:iDistance:round((subSx+1)/2)+iDistance,...
    round((subSy+1)/2)-iDistance:iDistance:round((subSy+1)/2)+iDistance) = 1;
subImSample1 = subImSample;

% horizontal lines: 8 10 14 20 32 48 64
subImSample = zeros(subSx, subSy);
iSpace = [4 6 8 10 14 20];
numLines = length(iSpace)+1;
ia = ones(1, numLines);
for i = 1:numLines-1
    ia(i+1) = ia(i) + iSpace(i);
end
subImSample(ia,:) = 1;
subImSample2 = subImSample;

% vertical lines
subImSample = zeros(subSx, subSy);
subImSample(:,ia) = 1;
subImSample3 = subImSample;

imOx = round((Sx-subSx*3)/2);
imOy = round((Sx-subSy)/2);
imSample(imOx:imOx+subSx-1,imOy:imOy+subSy-1) = subImSample1;
imSample(imOx+subSx:imOx+subSx*2-1,imOy:imOy+subSy-1) = subImSample2;
imSample(imOx+subSx*2:imOx+subSx*3-1,imOy:imOy+subSy-1) = flip(subImSample3,1);
WriteTifStack(imSample', 'imSample.tif', 32);
