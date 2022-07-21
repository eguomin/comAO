function [imgDs, zernPolynomials, Rc, waveFront_deltas] = zernretrieve_pre(imgs, r0, theta0, idx, p, c_delta, nFlag)
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
%   r0,theta0,idx: define the pupil aperture of the wavefront
%   1)r0: a 2D matrix of numbers between 0 and 1
%   2)theta0: a 2D matrix of angles (rad), has same size with r0
%   3)idx:
%       elements should be positive integers(>=1)
%   p: a vector of single indexes(Fringe convention) for Zernike components,
%   c_delta: a vector or matrix (known aberration) of Zernike coefficients
%   corresponding to p (should be normalized to phase unit: pi); if matrix,
%   the second dimension should be N-1

% By: Min Guo
% Jan 29, 2020
% Modifications: July 1, 2020
%   modify input pupil parameters(r, theta) from 1D vectors to 2D matrix 

r = r0(idx); % convert 2D matrix to vector with selected elements
theta = theta0(idx); % convert 2D matrix to vector with selected elements

[Sx, Sy, imgNum] = size(imgs); % image numbers:
deltaNum = size(c_delta,1);
if (imgNum - deltaNum)~= 1
    error('zernretrieve_pre:NMlength','deltaNum should be: imgNum -1.')
end
idxNum = length(r);
opNum = length(p);

imgDs = zeros(Sx, Sy, imgNum, 'single'); % FT of imgs
for k = 1:imgNum
    img = squeeze(imgs(:,:,k));
    imgD = fft2(img);
    imgDs(:,:,k) = imgD;
end    

zernPolynomials = create_zernpolybasis(p,r,theta,nFlag); % Zernike polynomial basis

% Calculate Rc
zernP2D = zeros(Sx, Sy,'single');
zernPhs = zeros(idxNum, opNum, 'single');
zernPvs = zeros(idxNum, opNum, 'single');
Rc = zeros(opNum, opNum, 'single');
zernP2Dtemp1 = zeros(Sx, Sy, 'single');
zernP2Dtemp2 = zeros(Sx, Sy, 'single');
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

waveFront_deltas = zeros(Sx,Sy,imgNum, 'single'); % known aberration, e.g., defocus
waveFront_delta = zeros(Sx,Sy, 'single');
for k = 2:imgNum
    waveFront_delta(idx) = create_wavefront(p, c_delta(k-1,:), r, theta, nFlag);
    waveFront_deltas(:,:,k) = waveFront_delta;
end

