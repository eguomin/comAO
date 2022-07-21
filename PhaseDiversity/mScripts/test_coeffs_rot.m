% A simple example: test coeffs_rot

% By: Min Guo, 
% Last update: Sep. 1, 2020

pIn = 1:32;
coeffs = zeros(1,32);
coeffs(3) = 0.3;

coeffsOut = coeffs_rot(coeffs, pIn, 45);
