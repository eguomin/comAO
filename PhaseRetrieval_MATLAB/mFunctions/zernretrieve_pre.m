function [imgDs, zernPolynomials, Rc, waveFront_deltas] = zernretrieve_pre(imgs, r, theta, idx, p, c_delta)
% perform pre-calculations for phase (zernike coefficients) retrieval
% update formula(recurrence): 
% Output
%   imgDs: FT of phase diversity images, third dimension: N, the number of images
%   zernPolynomials: basis of Zernike components
%   Rc: Rc matrix
%   waveFront_deltas: wavefronts corresponding to phase diversity images
%   (N)
% Input
%   imgs: phase diversity images, third dimension: N, the number of images
%   r,theta,idx: define the pupil aperture of the wavefront
%   1)r: a vector of numbers between 0 and 1
%   2)theta: a vector of angles, has same length with r
%   3)idx:
%       elements should be positive integers(>=1)
%   p: a vector of single indexes(Fringe convention) for Zernike components,
%   c_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (should be normalized to phase unit: pi); if matrix,
%   the second dimension should be N-1
% By: Min Guo
% Jan 29, 2020

[Sx, Sy, imgNum] = size(imgs); % image numbers: 2
deltaNum = size(c_delta,1);
if (imgNum - deltaNum)~= 1
    error('zernretrieve_pre:NMlength','deltaNum should be: imgNum -1.')
end
idxNum = length(r);
opNum = length(p);

imgDs = zeros(Sx, Sy, imgNum); % FT of imgs
for k = 1:imgNum
    img = squeeze(imgs(:,:,k));
    imgD = fft2(img);
    imgDs(:,:,k) = imgD;
end    

zernPolynomials = create_zernpolybasis(p,r,theta); % Zernike polynomial basis

% Calculate Rc
zernP2D = zeros(Sx,Sy);
zernPhs = zeros(idxNum,opNum);
zernPvs = zeros(idxNum,opNum);
Rc = zeros(opNum, opNum);
zernP2Dtemp1 = zeros(Sx,Sy);
zernP2Dtemp2 = zeros(Sx,Sy);
for m = 1:opNum
    zernP2D(idx) = zernPolynomials(:,m);
    zernP2Dtemp1(:,2:Sy) = zernP2D(:,1:Sy-1);
    zernP2Dshifted = zernP2Dtemp1 - zernP2D;
    zernPhs(:,m) = zernP2Dshifted(idx);
    zernP2Dtemp2(2:Sx,:) = zernP2D(1:Sx-1,:);
    zernP2Dshifted = zernP2Dtemp2 - zernP2D;
    zernPvs(:,m) = zernP2Dshifted(idx);
end
for m = 1:opNum
    for n = 1:opNum
        Rc(m,n) = dot(zernPhs(:,m),zernPhs(:,n)) + dot(zernPvs(:,m),zernPvs(:,n));
    end
end

waveFront_deltas = zeros(Sx,Sy,imgNum); % known aberration, e.g., defocus
waveFront_delta = zeros(Sx,Sy);
for k = 2:imgNum
    waveFront_delta(idx) = create_wavefront(p, c_delta(k-1,:), r, theta);
    waveFront_deltas(:,:,k) = waveFront_delta;
end

