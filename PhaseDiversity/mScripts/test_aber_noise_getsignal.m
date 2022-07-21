% test get signal
fileImgSample = ['..\..\..\computAO\Data\bootstrapping_new\none_astig\',...
    'Image_aberrated.tif'];
img0 = single(ReadTifStack(fileImgSample));
sType = 3;
threPerc = 0.1;
S = addnoise_getSignal(img0,sType, threPerc);
S
S/255*8192