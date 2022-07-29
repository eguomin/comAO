Control Software
===
---
LabVIEW (2018 32-bit) control software for the AO microscope, currently tested on the widefield setup:

1) master VIs for EMCCD control (widefield imaging);

2) master VIs for DM, WFS and EMCCD control, including data acquisition for DM calibration ;

3) master VIs for sensorless AO correction;

4) mini VIs for sample stage, laser shutter, DM and EMCCD etc;

5) sub VIs of sub-functions.


### Key VIs

- `AO_master_Widefield_v0.vi`: example to control EMCCD and perform widefield imaging (2D slice and z-stack).

- `AO_master_DM_WFS_v3.vi`: example to control DM, WFS and EMCCD.  Capable to:

    - set DM and do live wavefront measurement based on WFS;

    - perform widefield imaging;

    - check response of DM based on input of: a) zernike efficients and b) actuator voltages.

    - acquire data for DM calibration (to obtain command matrix) based on WFS or phase diversity.

    Note the bugs in this VI. 

- `AO_master_CorrLoop_PD_2D_v4.vi`: example to perform sensorless wavefront correction. This VI relies on the MATLAB function `recon_zern_fileIO` for wavefront reconstruction based on phase diversity.

- `AO_mini_stageLaserControl.vi`: mini example to control ASI stage, flip mirror and lasers.

- `AO_mini_AndorEMCCD.vi`: mini example to control EMCCD.

- `AO_mini_DM setup.vi`: mini example to set up DM.

- `AO_mini_wavefrontMeasure.vi`: mini example to measure wavefront.


To fully  run the VIs, it may require the installation of the drivers/SDK for the AO devices (Imagine Optics), EMCCD (Andor iXon), ASI stage (MS-2000), Thorlabs shutter, *etc*.