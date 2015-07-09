// matlab include
#include "mex.h"
#include "matrix.h"

// c++ include
#include "DT.h"
#include "SPBP.h"

#include <ctime>
#include <iostream>
#include <cmath>

#include "omp.h"
#include "openmpflag.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
// rhs (input): 
// rhs[0]: unary (n_state by n_node, single type)
// rhs[1]: params 
// lhs (output):
// lhs[0]: optimal state (n_node by 1)
{
	if ( nrhs != 2 ) { 
		mexPrintf("the number of input arg must be 2\n");
		return;
	}
	if ( ! mxIsSingle(prhs[0]) ) {
		mexPrintf("wrong unary potential format: must be single\n");
		return;
	}
	
	mwSize n_dim1 = mxGetNumberOfDimensions(prhs[0]);
	if ( n_dim1 != 2) {
		mexPrintf("wrong feature format -- it must be n_statexn_node\n");
        return;
	}
    size_t n_state, n_node;
    const mwSize* dims1 = mxGetDimensions(prhs[0]);
    n_state = dims1[0];
    n_node = dims1[1];
    
	// get parameter values
	mxArray *tmp;
	tmp = mxGetField(prhs[1], 0, "deformation_coeff");
	double deformation_coeff = mxGetScalar(tmp);

	tmp = mxGetField(prhs[1], 0, "truncate_const");
	double truncate_const = mxGetScalar(tmp);
	
	tmp = mxGetField(prhs[1], 0, "max_iter");
	int max_iter = (int)mxGetScalar(tmp);

	tmp = mxGetField(prhs[1], 0, "n_tree_level");
	int n_tree_level = (int)mxGetScalar(tmp);

	tmp = mxGetField(prhs[1], 0, "n_xstate");
	int n_xstate = (int)mxGetScalar(tmp);

	tmp = mxGetField(prhs[1], 0, "n_ystate");
	int n_ystate = (int)mxGetScalar(tmp);
	
	// run BP	
	SPBPParams params;
	params.n_sp_level = n_tree_level;
	params.n_state = n_state;
	params.n_xstate = n_xstate;
	params.n_ystate = n_ystate;
	SPBP bp(params);
	float* unary = (float*)mxGetData(prhs[0]);
	bp.RunBP(unary, deformation_coeff, truncate_const, max_iter);
	
	// output
	mwSize m = n_node;
	mwSize n = 1;
	plhs[0] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
	int* optimal_state = (int*)mxGetData(plhs[0]);
	bp.GetOptimalState(optimal_state);	
	
}
