                    Deformable Spatial Pyramid Matching for Fast Dense Correspondences

Install
This code is tested in Ubuntu 12.04 (GCC 4.4, 64bit) and Windows (VS2008 and 2012, 32bit) under MATLAB R2011a and R2013a.
To install, simply run mex_compile.m in MATLAB

Run example codes
- Before you run the codes, you have to set up the MATLAB path for the external library: run('vlfeat-0.9.17/toolbox/vl_setup.m')
DSPDemo.m : Include SIFT extraction (128-dim) and DSP matching 
DSPDemo_CT101.m: Include PCA-SIFT (20-dim) and DSP matching in Caltech-101 dataset
DSPDemo_LMO.m: Include PCA-SIFT (20-dim) and DSP matching in LabelMe Outdoor dataset

Multi-thread running
The codes run in multi-thread by default. 
To turn off the multi-thread option, open "openmpflag.h" and comment #define __OPEN_MP

Run-time
The run-time in the paper is from Intel Xeon X5690 (3.47GHz) with 16 cores.

Citation
Jaechul Kim, Ce Liu, Fei Sha, Kristen Grauman, 
Deformable Spatial Pyramid Matching for Fast Dense Correspondences, CVPR 2013







