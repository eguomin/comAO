tic;
lambda = 550; % nm
NA = 0.8;
Sx = 768;
Sy = Sx;
pixelSize = 162.5; % nm
[r, theta, idx] = def_pupilcoor(Sx, pixelSize, lambda, NA);
cTime1 = toc

p = 2:36;
zernPolynomials = create_zernpolybasis(p,r(idx),theta(idx));
cTime2 = toc


% Calculate Rc
opNum = length(p);
idxNum = length(r(idx));
zernP2D = zeros(Sx,Sy);
zernPhs = zeros(idxNum,opNum);
zernPvs = zeros(idxNum,opNum);
Rc = zeros(opNum, opNum);
for m = 1:opNum
    zernP2D(idx) = zernPolynomials(:,m);
    zernP2Dshifted = imtranslate(zernP2D,[1 0]) - zernP2D;
    zernPhs(:,m) = zernP2Dshifted(idx);
    zernP2Dshifted = imtranslate(zernP2D,[0 1]) - zernP2D;
    zernPvs(:,m) = zernP2Dshifted(idx);
end
cTime3 = toc
opNum = length(p);
idxNum = length(r(idx));
zernP2D = zeros(Sx,Sy);
zernPhs2 = zeros(idxNum,opNum);
zernPvs2 = zeros(idxNum,opNum);
Rc = zeros(opNum, opNum);
zernP2Dtemp1 = zeros(Sx,Sy);
zernP2Dtemp2 = zeros(Sx,Sy);
for m = 1:opNum
    zernP2D(idx) = zernPolynomials(:,m);
    zernP2Dtemp1(:,2:Sy) = zernP2D(:,1:Sy-1);
    zernP2Dshifted = zernP2Dtemp1 - zernP2D;
    zernPhs2(:,m) = zernP2Dshifted(idx);
    zernP2Dtemp2(2:Sx,:) = zernP2D(1:Sx-1,:);
    zernP2Dshifted = zernP2Dtemp2 - zernP2D;
    zernPvs2(:,m) = zernP2Dshifted(idx);
end
cTime4 = toc
cTime2 - cTime1
cTime3 - cTime2
cTime4 - cTime3
aa = zernPhs2 - zernPhs;
bb = zernPvs2 - zernPvs;
max(aa(:))
max(bb(:))