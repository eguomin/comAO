# computationalAO
Adaptive optics for fluorescence microscopy based on phase retrieval algorithm

Scripts and functions:
1) recon_zern.m: key function for estimating the zernike coefficients from phase diversity images and known phases (zernike coefficients). It calls two subfunctions -- zernretrieve_pre and zernretrieve_loop; 
2) test_recon_zern_simu.m: example to reconstruct aberrated phase (zernike coefficients) based on simulated images; 
2) test_recon_zern_exp.m: example to reconstruct aberrated phase (zernike coefficients) based on experimental images.
