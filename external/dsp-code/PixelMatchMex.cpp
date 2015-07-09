// matlab include
#include "mex.h"
#include "matrix.h"

// c++ include
#include "PixelMatch.h"

#include <ctime>
#include <cmath>
#include <string.h>

#include "omp.h"
#include "openmpflag.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
// rhs (input): 
// rhs[0]: feat1 (im1_height by im1_width by feat_dim, single type)
// rhs[1]: feat2 (im2_height by im2_width by feat_dim, single type)
// rhs[2]: pixel_disparity (2 by im1_height*im1_width -- order: height by width, int32 type) 
// rhs[3]: params
// lhs (output):
// lhs[0]: pixel_disparity
// lhs[1]: match_cost
{
	if ( nrhs != 4 ) { 
		mexPrintf("the number of input arg must be 4\n");
		return;
	}
	if ( ! mxIsSingle(prhs[0]) ||  ! mxIsSingle(prhs[1]) ) {
		mexPrintf("wrong descriptor format: must be single\n");
		return;
	}
	if ( ! mxIsInt32(prhs[2]) ) {
		mexPrintf("wrong disparity format: must be int32\n");
	}

	mwSize n_dim1 = mxGetNumberOfDimensions(prhs[0]);
	if ( n_dim1 != 3) {
		mexPrintf("wrong feature format -- it must be heightxwidthxfeat_dim\n");
        return;
	}
    size_t feat_dim1, height1, width1;
    const mwSize* dims1 = mxGetDimensions(prhs[0]);
    height1 = dims1[0];
    width1 = dims1[1];
    feat_dim1 = dims1[2];
	mwSize n_dim2 = mxGetNumberOfDimensions(prhs[1]);
	if ( n_dim2 != 3) {
		mexPrintf("wrong feature format -- it must be heightxwidthxfeat_dim\n");
        return;
	}
	size_t feat_dim2, height2, width2;
    const mwSize* dims2 = mxGetDimensions(prhs[1]);
    height2 = dims2[0];
    width2 = dims2[1];
    feat_dim2 = dims2[2];

	if ( feat_dim1 != feat_dim2 ){
		mexPrintf("different feature dimension\n");
		return;
	}
	size_t feat_dim = feat_dim1;
	
	// get descriptors
	float* feat1 = (float*)mxGetData(prhs[0]);
	float* feat2 = (float*)mxGetData(prhs[1]);

	// get pixel disparity
	int* pixel_disparity = (int*)mxGetData(prhs[2]);
	
	// get parameter values
	mxArray *tmp;
	tmp = mxGetField(prhs[3], 0, "search_grid_size");
	int search_grid_size = (int)mxGetScalar(tmp);

	tmp = mxGetField(prhs[3], 0, "search_radius");
	double* search_radius = mxGetPr(tmp);
	int search_radius_x = (int)search_radius[0];
	int search_radius_y = (int)search_radius[1];
	
	tmp = mxGetField(prhs[3], 0, "deform_coeff");
	double deform_coeff = mxGetScalar(tmp);

	tmp = mxGetField(prhs[3], 0, "truncate_const");
	double truncate_const = (float)mxGetScalar(tmp);
	
	// set parameters
	PixelMatchParams params;
	params.search_grid_size = search_grid_size;
	params.search_size_x = search_radius_x;
	params.search_size_y = search_radius_y;
	params.deformation_coeff = deform_coeff;
	params.truncate_const = truncate_const;

	// run the code
	PixelMatch pm;
	Descriptor descriptor1;
	descriptor1.desc = feat1;
	descriptor1.width = width1;
	descriptor1.height = height1;
	descriptor1.desc_dim = feat_dim;

	Descriptor descriptor2;
	descriptor2.desc = feat2;
	descriptor2.width = width2;
	descriptor2.height = height2;
	descriptor2.desc_dim = feat_dim;

	// output -- pixel's best disparity
	const mwSize* dims3 = mxGetDimensions(prhs[2]);
	plhs[0] = mxCreateNumericMatrix(dims3[0], dims3[1], mxINT32_CLASS, mxREAL);
	int* out_pixel_disparity = (int*)mxGetData(plhs[0]);
	memcpy(out_pixel_disparity, pixel_disparity, sizeof(int)*dims3[0]*dims3[1]);

	plhs[1] = mxCreateNumericMatrix(1, dims3[1], mxSINGLE_CLASS, mxREAL);
	float* match_cost = (float*)mxGetData(plhs[1]);

	pm.RunPixelMatch(descriptor1, descriptor2, out_pixel_disparity, match_cost, params);
	
}
