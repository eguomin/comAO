% test ipm

fileImgName = 'C:\Programs\computAO\Retrieval\MetricTest\FixedCell_groundTruth_crop.tif';
img0 = double(ReadTifStack(fileImgName));
img = max(img0-335,0);
mType = 1;
M1 = calculate_iqm(img,mType);
mType = 2;
M2 = calculate_iqm(img,mType);
mType = 3;
M3 = calculate_iqm(img,mType);
mType = 4;
M4 = calculate_iqm(img,mType);
mType = 6;
M5 = calculate_iqm(img,mType);
mType = 5;
M6 = calculate_iqm(img,mType);
mType = 7;
M7 = calculate_iqm(img,mType);


