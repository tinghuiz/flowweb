
#include "omp.h"
#include "openmpflag.h"
#include <algorithm>

#include "DT.h"


using namespace std;

/* dt of 1d function using squared distance */
void DT1d(float *f, int n, float c, float d) 
// v(x) = min(c|x|, d)
{
 	float min_val = *min_element(f, f+n);
	min_val += d;

	for (int i=1; i < n; i++) {
		f[i] = min(f[i], f[i-1]+c);
	}
	for (int i=n-2; i >=0 ; i--) {
		f[i] = min(f[i], f[i+1]+c);
	}
	for (int i=0; i < n; i++) {
		f[i] = min(f[i], min_val);
	}
}

void DT2d(float* im, int width, int height, float c, float d)
{
	
	// transform along columns
//  never do this!
//  #ifdef __OPEN_MP
// 		#pragma omp parallel for	
//  #endif
	//float *f = new float[max(height,width)];

	for (int x = 0; x < width; x++) {
		float *f = new float[height];
		for (int y = 0; y < height; y++) {
			int idx = y + x*height;
			f[y] = im[idx];
		}
		//int offset = x*height;
		//memcpy(im + offset, f, sizeof(float)*height);

		DT1d(f, height, c, d);

		for (int y = 0; y < height; y++) {
			int idx = y + x*height;
			im[idx] = f[y];
		}
		//memcpy(im + offset, f, sizeof(float)*height);
		delete [] f;
	}
    
//  never do this!	
// 	#ifdef __OPEN_MP
// 		#pragma omp parallel for	
//  #endif
	// transform along rows
	//#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		float *f = new float[width];
		for (int x = 0; x < width; x++) {
			int idx = y + x*height;
			f[x] = im[idx];
		}

		DT1d(f, width, c, d);
		
		for (int x = 0; x < width; x++) {
			int idx = y + x*height;		  
			im[idx] = f[x];
		}
		delete [] f;
	}
}

void DT3d(float* im, int width, int height, int depth, float c[3], float d[3])
// c[3], d[3] --> order of (height, width, depth)
{

	// xy-plane
//#ifdef __OPEN_MP
//	#pragma omp parallel for	
//#endif

	for ( int z =0 ; z < depth ; z++) {
		
		float* im_z = im + z*height*width;

		// transform along columns
		for (int x = 0; x < width; x++) {
			float *f = new float[height];
			for (int y = 0; y < height; y++) {
				int idx = y + x*height;
				f[y] = im_z[idx];
			}
			
			DT1d(f, height, c[0], d[0]);

			for (int y = 0; y < height; y++) {
				int idx = y + x*height;
				im_z[idx] = f[y];
			}
			delete [] f;
		}
	    
		// transform along rows
		for (int y = 0; y < height; y++) {
			float *f = new float[width];
			for (int x = 0; x < width; x++) {
				int idx = y + x*height;
				f[x] = im_z[idx];
			}

			DT1d(f, width, c[1], d[1]);
			
			for (int x = 0; x < width; x++) {
				int idx = y + x*height;		  
				im_z[idx] = f[x];
			}
			delete [] f;
		}
	}

	//along the z-axis
	int n_xy_elem = width*height;

//#ifdef __OPEN_MP
//	#pragma omp parallel for	
//#endif

	for ( int pos = 0 ; pos < n_xy_elem ; pos++) {

		float *f = new float[depth];
			
		for (int z = 0; z < depth; z++) {
			int idx = pos + z*n_xy_elem;
			f[z] = im[idx];
		}

		DT1d(f, depth, c[2], d[2]);
			
		for (int z = 0; z < depth; z++) {
			int idx = pos + z*n_xy_elem;
			im[idx] = f[z];
		}
		delete [] f;
	}
		
}

void Integral2D(float* src, float* dst, int width, int height)
{
	// first col only
	float rs = 0.0f;
	for(int i=0; i<height; i++) {
		rs += src[i]; 
		dst[i] = rs;
	}
	
	// remaining cells are sum above and to the left
	for(int i=1; i<width; i++) {

		rs = 0.0f;
		for(int j=0; j<height; j++) {
			rs += src[i*height+j]; 
			dst[i*height+j] = rs + dst[(i-1)*height+j];
		}
	}
}

void Integral2D(double* src, double* dst, int width, int height)
{
	// first col only
	double rs = 0.0;
	for(int i=0; i<height; i++) {
		rs += src[i]; 
		dst[i] = rs;
	}
	
	// remaining cells are sum above and to the left
	for(int i=1; i<width; i++) {

		rs = 0.0;
		for(int j=0; j<height; j++) {
			rs += src[i*height+j]; 
			dst[i*height+j] = rs + dst[(i-1)*height+j];
		}
	}
}