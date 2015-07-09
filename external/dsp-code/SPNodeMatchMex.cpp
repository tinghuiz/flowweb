// matlab include
#include "mex.h"
#include "matrix.h"

// c++ include
#include "DT.h"
#include "SPNodeMatch.h"

#include <ctime>
#include <iostream>
#include <cmath>

#include "omp.h"
#include "openmpflag.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
// rhs (input): 
// rhs[0]: feat1 (im1_height by im1_width by feat_dim, single type)
// rhs[1]: feat2 (im2_height by im2_width by feat_dim, single type)
// rhs[2]: params
// lhs (output):
// lhs[0]: node_cost (n_disp by n_node )
// lhs[1]: node disparity (n_disp by 1)
{
	if ( nrhs != 3 ) { 
		mexPrintf("the number of input arg must be 3\n");
		return;
	}
	if ( ! mxIsSingle(prhs[0]) ||  ! mxIsSingle(prhs[1]) ) {
		mexPrintf("wrong descriptor format: must be single\n");
		return;
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
	
	// get parameter values
	mxArray *tmp;
	tmp = mxGetField(prhs[2], 0, "n_tree_level");
	int n_tree_level = (int)mxGetScalar(tmp);

	tmp = mxGetField(prhs[2], 0, "search_radius");
	double* search_radius = mxGetPr(tmp);
	int search_radius_x = (int)search_radius[0];
	int search_radius_y = (int)search_radius[1];
	
	tmp = mxGetField(prhs[2], 0, "grid_size");
	int grid_size = (int)mxGetScalar(tmp);
	
	tmp = mxGetField(prhs[2], 0, "truncate_const");
	double truncate_const = mxGetScalar(tmp);
	
	// set parameters
	MatchParams params;
	params.sample_grid_size = grid_size;
	params.n_tree_level = n_tree_level;
	params.search_size_x = search_radius_x;
	params.search_size_y = search_radius_y;
	params.truncate_const = truncate_const;

	SPNodeMatch sp_match;
	int n_node = (int)pow(4.0,n_tree_level)/3;
	
	int dist_mat_width(0), dist_mat_height(0);
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int disp_x_range[2], disp_y_range[2];
	int n_disp_x(0), n_disp_y(0);
	
	ImSize im_size1, im_size2;
	im_size1.width = width1;
	im_size1.height = height1;
	im_size2.width = width2;
	im_size2.height = height2;

	sp_match.ComputeDistMatSize(im_size1, im_size2,
						grid_size,
						params.search_size_x, params.search_size_y,
						dist_mat_width, dist_mat_height,
						dist_mat_dx_range, dist_mat_dy_range,
						disp_x_range, disp_y_range,
						n_disp_x, n_disp_y);
	
	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;

	
	// output1 -- node's match cost
	mwSize m = n_disp;
	mwSize n = n_node;
	plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
	float* node_cost = (float*)mxGetData(plhs[0]);
	
	// output2 -- node's disparity
	m = 2;
	n = n_disp;
	plhs[1] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
	int* node_disparity = (int*)mxGetData(plhs[1]);

	// output3 -- auxilary info
	m = 1;
	n = 2;
	plhs[2] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
	int *aux_info = (int*)mxGetData(plhs[2]);
	aux_info[0] = n_dist_mat_dx;
	aux_info[1] = n_dist_mat_dy;


	

	// run match

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

	sp_match.NodeMatch(descriptor1, descriptor2, params, node_cost, node_disparity);
	
}
