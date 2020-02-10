function I = addpoissonnoise(I0)
% There should not be negetive values in I0
I = double(imnoise(uint16(I0), 'poisson' ));


