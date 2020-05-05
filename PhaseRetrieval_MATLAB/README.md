# Phase retrieval algorithm
MATLAB code for estimate of the zernike coefficients of aberration based on phase diversity images

Main scripts and functions:
1) recon_zern.m (in mFunctions folder): key function for estimating aberrated phase (Zernike coefficients) from phase diversity images and known phases (Zernike coefficients). It calls two subfunctions -- zernretrieve_pre and zernretrieve_loop; 
2) test_recon_zern_simu.m: example to reconstruct aberrated phase (Zernike coefficients) based on simulated images; 
3) test_recon_zern_exp.m: example to reconstruct aberrated phase (Zernike coefficients) based on experimental images.
4) test_recon_zern_fileIO.m: example to reconstruct aberrated phase based on experimental data directly from the microscope, by calling recon_zern_fileIO.m (in mFunctions folder). The recon_zern_fileIO.m is a key function that interprets the raw data/files to the reconstruction procedure by a series of prepocessing operations (background subtraction, image rotation/flipping, cropping and edge smoothing etc.).
