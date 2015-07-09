#pragma once

//#include "omp.h"
//#include "openmpflag.h"
//#include <string.h>
//#include "mex.h"

/* dt of 1d function using squared distance */
void DT1d(float *f, int n, float c, float d);
void DT2d(float* im, int width, int height, float c, float d);
void DT3d(float* im, int width, int height, int depth, float c[3], float d[3]);
void Integral2D(float* src, float* dst, int width, int height);
void Integral2D(double* src, double* dst, int width, int height);