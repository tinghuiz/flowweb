#include "PixelMatch.h"

#include <cfloat>
//#include <string.h>
//#include <functional>
//#include <numeric>
#include <algorithm>
//#include <iostream>
#include <cmath>

#include "omp.h"
#include "openmpflag.h"
#include "mex.h"

using namespace std;

void PixelMatch::RunPixelMatch(Descriptor descriptor1, 
							   Descriptor descriptor2,
							   int* bbox,
							   int* bbox_disparity,
							   int* pixel_disparity,
							   float* match_cost,
							   int n_bbox,
							   PixelMatchParams& params)
{

	// parameters
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int n_search(0);
	for ( int i=-params.search_size_x; i <= params.search_size_x; i += params.search_grid_size) {
		for ( int j=-params.search_size_y; j <= params.search_size_y;  j += params.search_grid_size) {
			n_search++;
		}
	}
	int* search_range = new int[2*n_search];
	int idx(0);
	for ( int i=-params.search_size_x; i <= params.search_size_x; i += params.search_grid_size) {
		for ( int j=-params.search_size_y; j <= params.search_size_y; j += params.search_grid_size) {
			search_range[2*idx] = i;
			search_range[2*idx+1] = j;
			idx++;
		}
	}
	
	float* deformation_cost = new float[n_search];
	for ( int i=0 ; i < n_search ; i++) {
		float dx = search_range[2*i];
		float dy = search_range[2*i+1];
		deformation_cost[i] = params.deformation_coeff*sqrt(dx*dx + dy*dy);
	}


	for ( int i=0 ; i < n_bbox; i++) {
		int* bbox_i = bbox + 4*i;
		int* ref_disparity_i = bbox_disparity + 2*i;

		UpdateDescDist(descriptor1, descriptor2, params,
		               match_cost, pixel_disparity,
					   bbox_i, ref_disparity_i,
					   search_range, deformation_cost, n_search);
	}
}

void PixelMatch::RunPixelMatch(Descriptor descriptor1, 
							   Descriptor descriptor2,
							   int* pixel_disparity, 
							   float* match_cost,
							   PixelMatchParams& params)
{
	// parameters
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int n_search(0);
	for ( int i=-params.search_size_x; i <= params.search_size_x; i += params.search_grid_size) {
		for ( int j=-params.search_size_y; j <= params.search_size_y;  j += params.search_grid_size) {
			n_search++;
		}
	}
	int* search_range = new int[2*n_search];
	int idx(0);
	for ( int i=-params.search_size_x; i <= params.search_size_x; i += params.search_grid_size) {
		for ( int j=-params.search_size_y; j <= params.search_size_y; j += params.search_grid_size) {
			search_range[2*idx] = i;
			search_range[2*idx+1] = j;
			idx++;
		}
	}
	
	// compute dist_mat
	//double start = omp_get_wtime();
	float* dist_mat = new float[width1*height1*n_search];
	int n_elem = width1*height1*n_search;
	for ( int i=0 ; i < n_elem ; i++) {
		dist_mat[i] = params.truncate_const;
	}
	ComputeDescDist(descriptor1, descriptor2, params,
					dist_mat, pixel_disparity, search_range, n_search);
	//double end = omp_get_wtime();
	//mexPrintf("%6.2f\n", end-start);

	// find the best disparity
	//start = omp_get_wtime();
	ComputeBestDisparity(dist_mat, pixel_disparity, match_cost,
						 search_range, params.deformation_coeff, 
						 width1, height1, n_search);
	//end = omp_get_wtime();
	//mexPrintf("%6.2f\n", end-start);


	delete [] dist_mat;

}

void PixelMatch::RunPixelMatch(Descriptor descriptor1, 
							   Descriptor descriptor2,
							   unsigned char* mask1,
							   unsigned char* mask2,
							   int* pixel_disparity, 
							   float* match_cost,
							   PixelMatchParams& params)
{
	// parameters
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int n_search(0);
	for ( int i=-params.search_size_x; i <= params.search_size_x; i += params.search_grid_size) {
		for ( int j=-params.search_size_y; j <= params.search_size_y;  j += params.search_grid_size) {
			n_search++;
		}
	}
	int* search_range = new int[2*n_search];
	int idx(0);
	for ( int i=-params.search_size_x; i <= params.search_size_x; i += params.search_grid_size) {
		for ( int j=-params.search_size_y; j <= params.search_size_y; j += params.search_grid_size) {
			search_range[2*idx] = i;
			search_range[2*idx+1] = j;
			idx++;
		}
	}
	
	// compute dist_mat
	//double start = omp_get_wtime();
	float* dist_mat = new float[width1*height1*n_search];
	int n_elem = width1*height1*n_search;
	for ( int i=0 ; i < n_elem ; i++) {
		dist_mat[i] = params.truncate_const;
	}
	ComputeDescDist(descriptor1, descriptor2, mask1, mask2, params,
					dist_mat, pixel_disparity, search_range, n_search);
	//double end = omp_get_wtime();
	//mexPrintf("%6.2f\n", end-start);

	// find the best disparity
	//start = omp_get_wtime();
	ComputeBestDisparity(dist_mat, pixel_disparity, match_cost,
						 search_range, params.deformation_coeff, 
						 width1, height1, n_search);
	//end = omp_get_wtime();
	//mexPrintf("%6.2f\n", end-start);


	delete [] dist_mat;

}


void PixelMatch::UpdateDescDist(Descriptor descriptor1,
								Descriptor descriptor2,
								PixelMatchParams& params,
								float* match_cost,
								int* pixel_disparity,
								int bbox[4],
								int ref_disparity[2],						
								int* search_range,
								float* deformation_cost,
								int n_search)
{
	float* desc1 = descriptor1.desc;
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int desc_dim = descriptor1.desc_dim;

	float* desc2 = descriptor2.desc;
	int width2 = descriptor2.width;
	int height2 = descriptor2.height;
	
	int n_pixel1 = width1*height1;
	int n_pixel2 = width2*height2;
	
	int ref_dx = ref_disparity[0];
	int ref_dy = ref_disparity[1];
	int lx = bbox[0];
	int rx = bbox[1];
	int ty = bbox[2];
	int dy = bbox[3];

#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_pixel1; pos++) {
		int x1 = pos/height1;
		int y1 = pos%height1;
		if ( x1 < lx || x1 > rx || y1 < ty || y1 > dy ) continue;

		int ptr1_desc_idx = pos;
		
		// for search range
		float min_dist(FLT_MAX);
		int min_idx(-1);
		for ( int i=0 ; i < n_search ; i++) {
			int dx = ref_dx + search_range[2*i];
			int dy = ref_dy + search_range[2*i+1];
		
			int x2 = x1 + dx;
			if ( x2 < 0 || x2 >= width2 ) continue;

			int y2 = y1 + dy;
			if ( y2 < 0 || y2 >= height2 ) continue;

			
			// compute descriptor distance
			float l2(0.0f);
			int ptr2_desc_idx = y2 + x2*height2;
			for  ( int dim=0 ; dim < desc_dim; dim++) {
				float val1 = desc1[ptr1_desc_idx + dim*n_pixel1];
				float val2 = desc2[ptr2_desc_idx + dim*n_pixel2];
				l2 += (val1-val2)*(val1-val2);
			}
			float dist_i = sqrt(l2) + deformation_cost[i];
			if ( min_dist > dist_i ) {
				min_dist = dist_i;
				min_idx = i;
			}
		}

		// update
		if ( match_cost[pos] > min_dist ) {
			match_cost[pos] = min_dist;
			pixel_disparity[2*pos] = ref_dx + search_range[2*min_idx];
			pixel_disparity[2*pos+1] = ref_dy + search_range[2*min_idx+1];
		}

	}
}


void PixelMatch::ComputeDescDist(Descriptor descriptor1, 
								 Descriptor descriptor2, 
								 PixelMatchParams& params,
								 float* dist_mat,
								 int* pixel_disparity,
								 int* search_range,
								 int n_search)
{
	
	float* desc1 = descriptor1.desc;
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int desc_dim = descriptor1.desc_dim;

	float* desc2 = descriptor2.desc;
	int width2 = descriptor2.width;
	int height2 = descriptor2.height;
	
	int n_pixel1 = width1*height1;
	int n_pixel2 = width2*height2;


#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_pixel1; pos++) {
		int x1 = pos/height1;
		int y1 = pos%height1;
		int ptr1_desc_idx = pos;

		int ref_dx = pixel_disparity[2*pos];
		int ref_dy = pixel_disparity[2*pos+1];

		// for search range
		for ( int i=0 ; i < n_search ; i++) {
			int dx = ref_dx + search_range[2*i];
			int dy = ref_dy + search_range[2*i+1];
		
			int x2 = x1 + dx;
			if ( x2 < 0 || x2 >= width2 ) continue;

			int y2 = y1 + dy;
			if ( y2 < 0 || y2 >= height2 ) continue;

			// dist_mat idx at (y1,x1,dy,dx) -- grid range
			int dist_mat_idx = pos + i*n_pixel1;
			
			// compute descriptor distance
			float l2(0.0f);
			int ptr2_desc_idx = y2 + x2*height2;
			for  ( int dim=0 ; dim < desc_dim; dim++) {
				float val1 = desc1[ptr1_desc_idx + dim*n_pixel1];
				float val2 = desc2[ptr2_desc_idx + dim*n_pixel2];
				l2 += (val1-val2)*(val1-val2);
			}
			dist_mat[dist_mat_idx] = sqrt(l2);

		}
	}
}

void PixelMatch::ComputeDescDist(Descriptor descriptor1, 
								 Descriptor descriptor2, 
								 unsigned char* mask1,
								 unsigned char* mask2,
								 PixelMatchParams& params,
								 float* dist_mat,
								 int* pixel_disparity,
								 int* search_range,
								 int n_search)
{
	
	float* desc1 = descriptor1.desc;
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int desc_dim = descriptor1.desc_dim;

	float* desc2 = descriptor2.desc;
	int width2 = descriptor2.width;
	int height2 = descriptor2.height;
	
	int n_pixel1 = width1*height1;
	int n_pixel2 = width2*height2;


#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_pixel1; pos++) {
		int x1 = pos/height1;
		int y1 = pos%height1;
		int ptr1_desc_idx = pos;

		int ref_dx = pixel_disparity[2*pos];
		int ref_dy = pixel_disparity[2*pos+1];

		// for search range
		for ( int i=0 ; i < n_search ; i++) {
			int dx = ref_dx + search_range[2*i];
			int dy = ref_dy + search_range[2*i+1];
		
			int x2 = x1 + dx;
			if ( x2 < 0 || x2 >= width2 ) continue;

			int y2 = y1 + dy;
			if ( y2 < 0 || y2 >= height2 ) continue;

			// dist_mat idx at (y1,x1,dy,dx) -- grid range
			int dist_mat_idx = pos + i*n_pixel1;
			int ptr2_desc_idx = y2 + x2*height2;

			if ( mask1[ptr1_desc_idx] == 2 || 
				 mask2[ptr2_desc_idx] == 2 || 
				 mask1[ptr1_desc_idx] == mask2[ptr2_desc_idx] ) {
				// compute descriptor distance
				float l2(0.0f);
				for  ( int dim=0 ; dim < desc_dim; dim++) {
					float val1 = desc1[ptr1_desc_idx + dim*n_pixel1];
					float val2 = desc2[ptr2_desc_idx + dim*n_pixel2];
					l2 += (val1-val2)*(val1-val2);
				}
				l2 = sqrt(l2);
				dist_mat[dist_mat_idx] = min(l2, params.truncate_const);
			}
			else {
				dist_mat[dist_mat_idx] = params.mask_mismatch_const;
			}

		}
	}
}

void PixelMatch::ComputeBestDisparity(float* dist_mat, 
									  int* pixel_disparity,
									  float* match_cost,
									  int* search_range,
									  float deformation_coeff,
							          int width, int height, int n_search)
{
	int n_pixel = width*height;
	float* deformation_cost = new float [n_search];
	for ( int i =0; i < n_search ; i++) {
		float dx = (float)search_range[2*i];
		float dy = (float)search_range[2*i+1];
		deformation_cost[i] = deformation_coeff*(fabs(dx) + fabs(dy));
	}


#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_pixel; pos++) {
		
		float min_val(FLT_MAX);
		int min_idx(-1);
		for ( int i=0 ; i < n_search ; i++) {

			int dist_mat_idx = pos + i*n_pixel;
			
			float val = dist_mat[dist_mat_idx] + deformation_cost[i];
			if ( val < min_val ) {
				min_val = val;
				min_idx = i;
			}
		}
		
		pixel_disparity[2*pos] += search_range[2*min_idx];
		pixel_disparity[2*pos+1] += search_range[2*min_idx+1];

		int idx = pos + min_idx*n_pixel;
		match_cost[pos] = dist_mat[idx];
	}
}