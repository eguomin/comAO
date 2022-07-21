function coeffsOut = fileIO_coeffs2lvtxt(fileTxtOut, p0, coeffs0, conType)
% convert the order list of Zernike polynomials (OSA indices) to  
% Wyant single-indices, and exported as txt file complied with LabVIEW
% Input
%   fileTxtOut: strings of the output path/file name;
%   p0: a vector of single-index zernike mode (Wyant convention)
%   coeffs0: a vector of Zernike coefficients defined by p0 
%   conType: the convention of single-index
%       'ANSI2ANSI': OSA/ANSI indices to OSA/ANSI indices; [default]
%       'ANSI2Wyant': OSA/ANSI indices to OSA/ANSI indices;
%       'Wyant2ANSI': Wyant indices to OSA/ANSI indices;
%       'Wyant2Wyant': Wyant indices to Wyant indices;
%   zernNumLimit: number of the output coefficients;

% By: Min Guo
% July 30, 2020
if(nargin==3)
    conType = 'ANSI2ANSI';
end
coeffsSigns = ones(1,50);
switch(conType)
    case 'ANSI2ANSI'
        zernNumLimit = 35;
        coeffsTemp = zernidx2idxcoeffs(p0, coeffs0, ... 
            'ANSI2Wyant');
        coeffs = zernidx2idxcoeffs(1:length(coeffsTemp), ...
            coeffsTemp, 'Wyant2ANSI', zernNumLimit);
        coeffsSigns([1,3, 6:7, 10:11, 15:17, 21:24, 28:31, 36:39, ...
            45:49]) = -1; 
    case 'ANSI2Wyant'
        zernNumLimit = 50;
        coeffs = zernidx2idxcoeffs(p0, coeffs0, 'ANSI2Wyant', zernNumLimit);
        coeffsSigns([2,5, 7, 10, 12, 14, 15, 17, 19, 21, 23, 26, 28, ... 
            30, 32, 34, 35, 37, 39, 41, 43, 45, 47, 50]) = -1; 
    case 'Wyant2ANSI'
        zernNumLimit = 35;
        coeffs = zernidx2idxcoeffs(p0, coeffs0, 'Wyant2ANSI', zernNumLimit);
        coeffsSigns([1,3, 6:7, 10:11, 15:17, 21:24, 28:31, 36:39, ...
            45:49]) = -1; 
    case 'Wyant2Wyant'
        zernNumLimit = 50;
        coeffsTemp = zernidx2idxcoeffs(p0, coeffs0, ... 
            'Wyant2ANSI');
        coeffs = zernidx2idxcoeffs(1:length(coeffsTemp), ...
            coeffsTemp, 'ANSI2Wyant', zernNumLimit);
        coeffsSigns([2,5, 7, 10, 12, 14, 15, 17, 19, 21, 23, 26, 28, ... 
            30, 32, 34, 35, 37, 39, 41, 43, 45, 47, 50]) = -1; 
    otherwise
        error('fileIO_coeffs2lvtxt: wrong convetion type');
end
coeffs = -coeffs;
coeffsOut = coeffs;
zernNumMin = min(length(coeffs), length(coeffsSigns));
coeffsOut(1:zernNumMin) = coeffs(1:zernNumMin).*coeffsSigns(1:zernNumMin);
dlmwrite(fileTxtOut, coeffsOut ,'delimiter','\t','precision','%.6f');