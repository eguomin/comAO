function zernPolynomials = create_zernpolybasis(p,r,theta,nflag)
% create_zernpolybasis is to creat Zernike polynomial basis based on 
% singel-index in Fringe convention

% Output
%   zernPolynomils: a vector of wavefront for every (r,theta) pair position
% Input
%   p: a vector of single indexes(Fringe convention) for Zernike components,
%       elements should be positive integers(>=1)
%   r: a vector of numbers between 0 and 1
%   theta: a vector of angles, has same length with r
%   nflag: optional, nflag = 'norm' is corresponding to the normalized Zernike
%       functions
% By: Min Guo
% July 26, 2017

[n, m] = zernfringe2nm(p);

switch nargin
    case 3
        zernPolynomials = zernfun(n,m,r,theta);
    case 4
        zernPolynomials = zernfun(n,m,r,theta,nflag);
    otherwise
        error('zernfun2:nargin','Incorrect number of inputs.')
end
    
