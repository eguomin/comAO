function waveFront = create_wavefront(p,coeffs,r,theta,nflag)
% CREATWAVEFRONT is to creat wavefront based on Zernike coefficients with
% singel-index in Fringe convention
% Output
%   waveFront: a vector of wavefront for every (r,theta) pair position
% Input
%   p: a vector of single indexes(Fringe convention) for Zernike components,
%       elements should be positive integers(>=1)
%   coeffs: a vector of Zernike coefficients corresponding to p
%   r: a vector of numbers between 0 and 1
%   theta: a vector of angles, has same length with r
%   nflag: optional, nflag = 'norm' is corresponding to the normalized Zernike
%       functions
% By: Min Guo
% Dec 09, 2016

if length(p)~=length(coeffs)
    error('creatwavefront:NMlength','p and coeffs must be the same length.')
end

[n, m] = zernfringe2nm(p);

switch nargin
    case 4
        z = zernfun(n,m,r,theta);
    case 5
        z = zernfun(n,m,r,theta,nflag);
    otherwise
        error('zernfun2:nargin','Incorrect number of inputs.')
end
waveFront = zeros(size(r));
for i = 1:length(p)
    waveFront = waveFront+coeffs(i)*z(:,i);
end
    
