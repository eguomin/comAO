function p = zernnm2idx(n,m, conType)
% p = zernnm2idx(n,m, conType) converts the order list for Zernike polynomials
% from dual-index(n, m) to single-index p.
% 
% Output
%   p: single-index scale or vector of integer numbers.
% Input
%   n: radial order, scale or vector of integer numbers.
%   m: angular frequency, scale or vector of integer numbers.
%   conType: the convention of single-index
%       'ANSI': OSA/ANSI single-index; [default]
%           Thibos, L. N.; Applegate, R. A.; Schwiegerling, J. T.; Webb, R. (2002).
%           "Standards for reporting the optical aberrations of eyes" (PDF). 
%           Journal of Refractive Surgery. 18 (5): S652–60
%       'Fringe': Fringe indices, 'Wyant' + 1;
%       'Wyant': Wyant indices , used in HASO;
% Note: a LUT can also be created to accelerate the converting.

% By: Min Guo
% July 22, 2020
if(nargin==2)
    conType = 'ANSI'; % defult as OSA/ANSI convention
end

switch conType
    case 'ANSI'
        p = (n.*(n+2) + m)/2;
    case 'Fringe'
        p = (1+(n+abs(m))/2).^2 - 2*abs(m) + (1-sign(m))/2;
        p = floor(p);
    case 'Wyant'
        p = (1+(n+abs(m))/2).^2 - 2*abs(m) + (1-sign(m))/2 - 1;
        p = floor(p);
    otherwise
        error('zernnm2idx: wrong convetion type');
end