Phase diversity algorithm
===
---

MATLAB code for estimating the zernike coefficients of aberration based on phase diversity images

### 1. Test scripts

- `test_wavefront_2_coeffs.m`: example to convert Zernike coefficients to wavefront, and convert wavefront back to Zernike coefficients

- `test_recon_simu.m`: example to reconstruct aberrated phase (Zernike coefficients) based on synthetic data.    

- `test_recon_exp.m`: example to reconstruct aberrated phase (Zernike coefficients) based on experimental data.

- `test_reconDM_exp_batch.m`: example to perform phase diversity algorithm on image stack (defocus diversity phase)  for deformable mirror calibration - batch processing. It requires images preprocessed using ImageJ script `PreProcess.ijm`.

- `test_calc_command_matrix.m`: build interaction matrix and calculate command matrix.

- `test_calibration_comparison.m`: Check and compare DM response based on zernike input for calculated command matrixes.


### 2. Key functions

- `mFunctions/recon_zern.m`: key function for estimating aberrated phase (Zernike coefficients) from phase diversity images and known phases based on Zernike coefficients. It calls two subfunctions -- zernretrieve_pre and zernretrieve_loop.

- `mFunctions/recon_zern_fileIO.m`: key function for estimating aberrated phase (Zernike coefficients) from experimental data. It calls recon_zern.m and a few functions that handle file IO, conversion and other preprocessing operations.

- `mFunctions/recon_zern_custom.m`: key function for estimating aberrated phase (Zernike coefficients) from phase diversity images and known phases - The known phases are direct wavefront input instead of Zernike coefficients.