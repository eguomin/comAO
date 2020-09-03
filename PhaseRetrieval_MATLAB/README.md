# Phase retrieval algorithm
MATLAB code for estimating the zernike coefficients of aberration based on phase diversity images

Main scripts and functions:
1) test_recon_simu.m: example to reconstruct aberrated phase (Zernike coefficients) based on synthetic data; 
2) test_recon_exp.m: example to reconstruct aberrated phase (Zernike coefficients) based on experimental data.
3) mFunction/recon_zern.m: key function for estimating aberrated phase (Zernike coefficients) from phase diversity images and known phases (Zernike coefficients). It calls two subfunctions -- zernretrieve_pre and zernretrieve_loop; 
4) mFunction/recon_zern_fileIO.m: key function for estimating aberrated phase (Zernike coefficients) from experimental data. It calls recon_zern.m and a few functions that handle file IO, conversion and other preprocessing operations.
