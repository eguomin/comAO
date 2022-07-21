function [coeffs, coeffsLv] = fileIO_lvtxt2idx(fileZernIn, conType, zernNumLimit)
% import Zernike polynomials from txt file(Wyant indices) to  
% Wyant or OSA/ANSI indices
% Input
%   fileTxtIn: strings of the input path/file name;
%   conType: the convention of single-index
%       'ANSI2ANSI': OSA/ANSI indices to OSA/ANSI indices; [default]
%       'ANSI2Wyant': OSA/ANSI indices to OSA/ANSI indices;
%       'Wyant2ANSI': Wyant indices to OSA/ANSI indices;
%       'Wyant2Wyant': Wyant indices to Wyant indices;
%   zernNumLimit: number of the output coefficients; [default 28 to Mirao limit]

% By: Min Guo
% July 30, 2020
if(nargin==1)
    conType = 'ANSI2ANSI';
    zernNumLimit = 28;
elseif(nargin==2)
    zernNumLimit = 28;
end
coeffsRaw = importdata(fileZernIn);
[vectorNum, zernNum] = size(coeffsRaw);
coeffsSigns = ones(1,50); % sin term = -1; 
                   % plus spherical = 1 or -1 (alternately by orders
switch(conType)
    case 'ANSI2ANSI'
        coeffsSigns([1,3, 6:7, 10:11, 15:17, 21:24, 28:31, 36:39, ...
            45:49]) = -1; 
    case 'ANSI2Wyant'
        coeffsSigns([1,3, 6:7, 10:11, 15:17, 21:24, 28:31, 36:39, ...
            45:49]) = -1; 
    case 'Wyant2ANSI'
        coeffsSigns([2,5, 7, 10, 12, 14, 15, 17, 19, 21, 23, 26, 28, ... 
            30, 32, 34, 35, 37, 39, 41, 43, 45, 47, 50]) = -1;  
    case 'Wyant2Wyant'
        coeffsSigns([2,5, 7, 10, 12, 14, 15, 17, 19, 21, 23, 26, 28, ... 
            30, 32, 34, 35, 37, 39, 41, 43, 45, 47, 50]) = -1;  
    otherwise
end
coeffsOut = coeffsRaw;
zernNumMin = min(zernNum, length(coeffsSigns));
for i = 1:vectorNum
    coeffsOut(i,1:zernNumMin) = coeffsRaw(i,1:zernNumMin) ...
        .*coeffsSigns(1,1:zernNumMin);
end

coeffs = zeros(vectorNum, zernNumLimit);

for i = 1:vectorNum
    coeffs0 = coeffsOut(i,:);
    switch(conType)
        case 'ANSI2ANSI'
            coeffsTemp = zernidx2idxcoeffs(1:zernNum, coeffs0, ... 
                'ANSI2Wyant');
            coeffs(i,:) = zernidx2idxcoeffs(1:length(coeffsTemp), ...
                coeffsTemp, 'Wyant2ANSI', zernNumLimit);
        case 'ANSI2Wyant'
            coeffs(i,:) = zernidx2idxcoeffs(1:zernNum, coeffs0, ... 
                'ANSI2Wyant', zernNumLimit);
        case 'Wyant2ANSI'
            coeffs(i,:) = zernidx2idxcoeffs(1:zernNum, coeffs0, ... 
                'Wyant2ANSI', zernNumLimit);
        case 'Wyant2Wyant' 
            coeffsTemp = zernidx2idxcoeffs(1:zernNum, coeffs0, ... 
                'Wyant2ANSI');
            coeffs(i,:) = zernidx2idxcoeffs(1:length(coeffsTemp), ...
                coeffsTemp, 'ANSI2Wyant', zernNumLimit);
        otherwise
    end
end
coeffs = -single(coeffs);
coeffsLv = coeffsRaw(:,1:zernNumLimit);