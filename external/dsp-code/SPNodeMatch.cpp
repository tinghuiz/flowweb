#include "SPNodeMatch.h"
#include "DT.h"
#include "MyUtil.h"
#include "math.h"

#include <cfloat>
#include <string.h>
#include <functional>
#include <numeric>
#include <algorithm>
#include <iostream>

#include "omp.h"
#include "openmpflag.h"
#include "mex.h"

void SPNodeMatch::ComputeDistMatSize(ImSize im_size1, ImSize im_size2,
								  int grid_size, 
								  int search_size_x, int search_size_y,
								  int& dist_mat_width, int& dist_mat_height,
								  int dist_mat_dx_range[2], int dist_mat_dy_range[2],
								  int disp_x_range[2], int disp_y_range[2],
								  int &n_disp_x, int &n_disp_y)
{
	int width1 = im_size1.width;
	int height1 = im_size1.height;
	int width2 = im_size2.width;
	int height2 = im_size2.height;

	// x,y dimension
	dist_mat_width = 0;
	for  ( int i = 0 ; i < width1 ;  i += grid_size) {
		dist_mat_width++;
	}
	dist_mat_height = 0;
	for  ( int i = 0 ; i < height1 ;  i += grid_size) {
		dist_mat_height++;
	}

	int max_dx = grid_size*(search_size_x/grid_size);
	int min_dx = -max_dx;

	/*int min_dx(0); 
	int max_x1 = grid_size*(width1/grid_size); 
	for ( int x2=0 ; x2 < width2; x2 += grid_size) {
		int dx = x2 - max_x1;	
		if ( dx > -search_size_x ) {
			min_dx = dx;
			break;
		}
	}*/
	int max_dy = grid_size*(search_size_y/grid_size);
	int min_dy = -max_dy;
	/*int min_dy(0); 
	int max_y1 = grid_size*(height1/grid_size); 
	for ( int y2=0 ; y2 < height2; y2 += grid_size) {
		int dy = y2 - max_y1;	
		if ( dy > -search_size_y ) {
			min_dy = dy;
			break;
		}
	}*/
	dist_mat_dx_range[0] = min_dx;
	dist_mat_dx_range[1] = max_dx;
	dist_mat_dy_range[0] = min_dy;
	dist_mat_dy_range[1] = max_dy;


	// disparity range
	disp_x_range[0] = -(width1-1);
	disp_x_range[1] = width2-1;
	disp_y_range[0] = -(height1-1);
	disp_y_range[1] = height2-1;

	// disparity dimension
	n_disp_x = 0;
	for  ( int i = disp_x_range[0] ; i <= disp_x_range[1] ;  i += grid_size) {
		n_disp_x++;
	}
	n_disp_y = 0;
	for  ( int i = disp_y_range[0] ; i <= disp_y_range[1] ;  i += grid_size) {
		n_disp_y++;
	}
	
}

void SPNodeMatch::ComputeDescDist(Descriptor descriptor1, Descriptor descriptor2, 
							    int grid_size,
							    MatchParams& params,
							    float* dist_mat,
								int* disparity)
// desc1: y,x,dim (height1 by width1 by desc_dim)
// desc2: y,x,dim (height2 by width2 by desc_dim)
// dist_mat: y,x,dy,dx (grid range, but only selected disparities)
{
	//double start = omp_get_wtime();
	
	float* desc1 = descriptor1.desc;
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int desc_dim = descriptor1.desc_dim;

	float* desc2 = descriptor2.desc;
	int width2 = descriptor2.width;
	int height2 = descriptor2.height;
	
	ImSize im_size1, im_size2;
	im_size1.height = height1;
	im_size1.width = width1;
	im_size2.height = height2;
	im_size2.width = width2;

	int n_pixel1 = width1*height1;
	int n_pixel2 = width2*height2;
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;

	int dist_mat_width(0), dist_mat_height(0);
	int disp_x_range[2], disp_y_range[2];
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int n_disp_x(0), n_disp_y(0);

	ComputeDistMatSize(im_size1, im_size2,
		               grid_size,
					   search_size_x, search_size_y,
					   dist_mat_width, dist_mat_height,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);

	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;
			
#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_dist_mat_xy_elem; pos++) {
		int x_idx = pos/dist_mat_height;
		int y_idx = pos%dist_mat_height;
		int x1 = x_idx*grid_size;
		int y1 = y_idx*grid_size;
		int ptr1_desc_idx = y1 + x1*height1;

		// for each disparity
		for (int i=0 ; i < n_disp ; i++) {
			
			int dx = disparity[2*i];
			int dy = disparity[2*i+1];

			int x2 = x1 + dx;
			if ( x2 < 0 || x2 >= width2 ) continue;

			int y2 = y1 + dy;
			if ( y2 < 0 || y2 >= height2 ) continue;

			// dist_mat idx at (y1,x1,dy,dx) -- grid range
			int dist_mat_idx = pos + i*n_dist_mat_xy_elem;
			
			// compute descriptor distance
			float l2(0.0f);
			int ptr2_desc_idx = y2 + x2*height2;
			for  ( int dim=0 ; dim < desc_dim; dim++) {
				float val1 = desc1[ptr1_desc_idx + dim*n_pixel1];
				float val2 = desc2[ptr2_desc_idx + dim*n_pixel2];
				l2 += (val1-val2)*(val1-val2);
			}
			l2 = sqrt(l2);
			dist_mat[dist_mat_idx] = min(l2, params.truncate_const);

		}
	}
}

void SPNodeMatch::ComputeDescDist(Descriptor descriptor1, Descriptor descriptor2, 
							    int grid_size, float scale,
							    MatchParams& params,
							    float* dist_mat,
								int* disparity)
// desc1: y,x,dim (height1 by width1 by desc_dim)
// desc2: y,x,dim (height2 by width2 by desc_dim)
// dist_mat: y,x,dy,dx (grid range, but only selected disparities)
{
	//double start = omp_get_wtime();
	
	float* desc1 = descriptor1.desc;
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int desc_dim = descriptor1.desc_dim;
	int cx1 = width1/2;
	int cy1 = height1/2;

	float* desc2 = descriptor2.desc;
	int width2 = descriptor2.width;
	int height2 = descriptor2.height;
	
	ImSize im_size1, im_size2;
	im_size1.height = height1;
	im_size1.width = width1;
	im_size2.height = height2;
	im_size2.width = width2;

	int n_pixel1 = width1*height1;
	int n_pixel2 = width2*height2;
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;

	int dist_mat_width(0), dist_mat_height(0);
	int disp_x_range[2], disp_y_range[2];
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int n_disp_x(0), n_disp_y(0);

	ComputeDistMatSize(im_size1, im_size2,
		               grid_size,
					   search_size_x, search_size_y,
					   dist_mat_width, dist_mat_height,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);

	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;
			
#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_dist_mat_xy_elem; pos++) {
		int x_idx = pos/dist_mat_height;
		int y_idx = pos%dist_mat_height;
		int x1 = x_idx*grid_size;
		int y1 = y_idx*grid_size;
		int ptr1_desc_idx = y1 + x1*height1;
		int dx_offset = scale*(x1-cx1) + cx1;
		int dy_offset = scale*(y1-cy1) + cy1;
		// for each disparity
		for (int i=0 ; i < n_disp ; i++) {
			
			int dx = disparity[2*i];
			int dy = disparity[2*i+1];

			int x2 = dx_offset + dx;
			if ( x2 < 0 || x2 >= width2 ) continue;

			int y2 = dy_offset + dy;
			if ( y2 < 0 || y2 >= height2 ) continue;

			// dist_mat idx at (y1,x1,dy,dx) -- grid range
			int dist_mat_idx = pos + i*n_dist_mat_xy_elem;
			
			// compute descriptor distance
			float l2(0.0f);
			int ptr2_desc_idx = y2 + x2*height2;
			for  ( int dim=0 ; dim < desc_dim; dim++) {
				float val1 = desc1[ptr1_desc_idx + dim*n_pixel1];
				float val2 = desc2[ptr2_desc_idx + dim*n_pixel2];
				l2 += (val1-val2)*(val1-val2);
			}
			l2 = sqrt(l2);
			dist_mat[dist_mat_idx] = min(l2, params.truncate_const);

		}
	}
}

void SPNodeMatch::NodeMatch(Descriptor descriptor1, 
							Descriptor descriptor2,
						    MatchParams& params, 
							float scale,
						    float* node_cost,
							int* node_disparity)
{
	ImSize im_size1;
	im_size1.width = descriptor1.width;
	im_size1.height = descriptor1.height;

	ImSize im_size2;
	im_size2.width = descriptor2.width;
	im_size2.height = descriptor2.height;

	int grid_size = params.sample_grid_size;
	int n_tree_level = params.n_tree_level;
	int n_node = ((int)pow(4.0,n_tree_level)-1)/3;
	
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;
		
	int dist_mat_width(0), dist_mat_height(0);
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int disp_x_range[2], disp_y_range[2];
	int n_disp_x(0), n_disp_y(0);
	
	ComputeDistMatSize(im_size1, im_size2,
					   grid_size,
					   search_size_x, search_size_y,
					   dist_mat_width, dist_mat_height,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);

	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;
	int disp_min_x = dist_mat_dx_range[0];
	int disp_min_y = dist_mat_dy_range[0];

	// allocate memory - dist_mat
	int n_elem = dist_mat_width*dist_mat_height*n_dist_mat_dx*n_dist_mat_dy;
	float* dist_mat = new float[n_elem];
	for ( int i=0 ; i < n_elem; i++) {
		dist_mat[i] = params.truncate_const;
	}
	// allocate memory - disparity
	for ( int dx=0; dx < n_dist_mat_dx ; dx++) {
		for ( int dy=0; dy < n_dist_mat_dy; dy++) {
			int x = dx*grid_size + disp_min_x;
			int y = dy*grid_size + disp_min_y;
			int idx = dy + dx*n_dist_mat_dy;
			node_disparity[2*idx] = x;
			node_disparity[2*idx+1] = y;
		}
	}

	
	ComputeDescDist(descriptor1, descriptor2, 
	    			grid_size, scale, params, dist_mat, node_disparity);


	// compute node's match cost	
	ComputeNodeCost(dist_mat,
		            im_size1, im_size2, 
					params,
					grid_size, 
					node_cost);

	delete [] dist_mat;	
}

void SPNodeMatch::NodeMatch(Descriptor descriptor1, 
							Descriptor descriptor2,
						    MatchParams& params, 
						    float* node_cost,
							int* node_disparity)
{
	ImSize im_size1;
	im_size1.width = descriptor1.width;
	im_size1.height = descriptor1.height;

	ImSize im_size2;
	im_size2.width = descriptor2.width;
	im_size2.height = descriptor2.height;

	int grid_size = params.sample_grid_size;
	int n_tree_level = params.n_tree_level;
	int n_node = ((int)pow(4.0,n_tree_level)-1)/3;
	
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;
		
	int dist_mat_width(0), dist_mat_height(0);
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int disp_x_range[2], disp_y_range[2];
	int n_disp_x(0), n_disp_y(0);
	
	ComputeDistMatSize(im_size1, im_size2,
					   grid_size,
					   search_size_x, search_size_y,
					   dist_mat_width, dist_mat_height,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);

	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;
	int disp_min_x = dist_mat_dx_range[0];
	int disp_min_y = dist_mat_dy_range[0];

	// allocate memory - dist_mat
	int n_elem = dist_mat_width*dist_mat_height*n_dist_mat_dx*n_dist_mat_dy;
	float* dist_mat = new float[n_elem];
	for ( int i=0 ; i < n_elem; i++) {
		dist_mat[i] = params.truncate_const;
	}
	// allocate memory - disparity
	for ( int dx=0; dx < n_dist_mat_dx ; dx++) {
		for ( int dy=0; dy < n_dist_mat_dy; dy++) {
			int x = dx*grid_size + disp_min_x;
			int y = dy*grid_size + disp_min_y;
			int idx = dy + dx*n_dist_mat_dy;
			node_disparity[2*idx] = x;
			node_disparity[2*idx+1] = y;
		}
	}

	
	ComputeDescDist(descriptor1, descriptor2, 
	    			grid_size, params, dist_mat, node_disparity);


	// compute node's match cost	
	ComputeNodeCost(dist_mat,
		            im_size1, im_size2, 
					params,
					grid_size, 
					node_cost);

	delete [] dist_mat;	
}

void SPNodeMatch::ComputeNodeCost(
					  float* dist_mat, 
					  ImSize im_size1, ImSize im_size2,
					  MatchParams& params, 
					  int grid_size,
					  float* node_cost)
// dist_mat: y,x,dy,dx (only for selected disparity)
// disp_mask: same as dist_mat
// disparity: dx,dy for selected disparity
{
	//double start = omp_get_wtime();

	int width1 = im_size1.width;
	int height1 = im_size1.height;
	int width2 = im_size2.width;
	int height2 = im_size2.height;

	int n_tree_level = params.n_tree_level;
	int n_node = ((int)pow(4.0,n_tree_level)-1)/3;
	
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;
	
	int dist_mat_w(0), dist_mat_h(0), n_disp_x(0), n_disp_y(0);
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int disp_x_range[2], disp_y_range[2];
	ComputeDistMatSize(im_size1, im_size2,
					   grid_size,
					   search_size_x, search_size_y,
					   dist_mat_w, dist_mat_h,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);
	
	int n_dist_mat_xy = dist_mat_w*dist_mat_h;
	int min_dx = dist_mat_dx_range[0];
	int min_dy = dist_mat_dy_range[0];
	
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;

	/* compute node cost */
	for (int disp=0; disp < n_disp; disp++) {
	
		// integral image 
		float* dist_integral = new float[n_dist_mat_xy];
		Integral2D(dist_mat + n_dist_mat_xy*disp, dist_integral, dist_mat_w, dist_mat_h);
		
		// for tree level
		for (int i=0; i < n_tree_level ; i++) { 		
		
			int node_idx_offset = ((int)pow(4.0,i)-1)/3;
			int n_bin_1d = (int)pow(2.0,i);
			int x_bin_size = (int)((double)dist_mat_w/(double)n_bin_1d);			
			int y_bin_size = (int)((double)dist_mat_h/(double)n_bin_1d);	

			// for each node
			for (int bx=0 ; bx < n_bin_1d; bx++) {
				int lx = bx*x_bin_size;
				int rx = (bx+1)*x_bin_size-1;
				for (int by=0 ; by < n_bin_1d ; by++) {
				
					int ty = by*y_bin_size;
					int dy = (by+1)*y_bin_size-1;
					
					int node_idx = node_idx_offset + by + bx*n_bin_1d;
					
					int idx1 = ty + lx*dist_mat_h;
					int idx2 = ty + rx*dist_mat_h;
					int idx3 = dy + lx*dist_mat_h;
					int idx4 = dy + rx*dist_mat_h;
					
					node_cost[node_idx*n_disp + disp] =
						(dist_integral[idx1] + dist_integral[idx4] - dist_integral[idx2] - dist_integral[idx3]);
					node_cost[node_idx*n_disp + disp] /= (dy-ty+1)*(rx-lx+1);

				}
			}
		}
		delete [] dist_integral;
	}

}


void SPNodeMatch::NodeMatch(Descriptor descriptor1, 
						    Descriptor descriptor2,
						    unsigned char* mask1,
						    unsigned char* mask2,
						    MatchParams& params,
						    float* node_cost,
						    int* node_disparity)
{
	ImSize im_size1;
	im_size1.width = descriptor1.width;
	im_size1.height = descriptor1.height;

	ImSize im_size2;
	im_size2.width = descriptor2.width;
	im_size2.height = descriptor2.height;

	int grid_size = params.sample_grid_size;
	int n_tree_level = params.n_tree_level;
	int n_node = ((int)pow(4.0,n_tree_level)-1)/3;
	
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;
		
	int dist_mat_width(0), dist_mat_height(0);
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int disp_x_range[2], disp_y_range[2];
	int n_disp_x(0), n_disp_y(0);
	
	ComputeDistMatSize(im_size1, im_size2,
					   grid_size,
					   search_size_x, search_size_y,
					   dist_mat_width, dist_mat_height,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);

	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;
	int disp_min_x = dist_mat_dx_range[0];
	int disp_min_y = dist_mat_dy_range[0];

	// allocate memory - dist_mat
	int n_elem = dist_mat_width*dist_mat_height*n_dist_mat_dx*n_dist_mat_dy;
	float* dist_mat = new float[n_elem];
	for ( int i=0 ; i < n_elem; i++) {
		dist_mat[i] = params.truncate_const;
	}
	// allocate memory - disparity
	for ( int dx=0; dx < n_dist_mat_dx ; dx++) {
		for ( int dy=0; dy < n_dist_mat_dy; dy++) {
			int x = dx*grid_size + disp_min_x;
			int y = dy*grid_size + disp_min_y;
			int idx = dy + dx*n_dist_mat_dy;
			node_disparity[2*idx] = x;
			node_disparity[2*idx+1] = y;
		}
	}

	
	ComputeDescDist(descriptor1, descriptor2, mask1, mask2,
	    			grid_size, params, dist_mat, node_disparity);


	// compute node's match cost	
	ComputeNodeCost(dist_mat,
		            im_size1, im_size2, 
					params,
					grid_size, 
					node_cost);

	delete [] dist_mat;	
}

void SPNodeMatch::ComputeDescDist(Descriptor descriptor1, Descriptor descriptor2, 
								  unsigned char* mask1, unsigned char* mask2,
								  int grid_size,
								  MatchParams& params,
								  float* dist_mat,
								  int* disparity)
// desc1: y,x,dim (height1 by width1 by desc_dim)
// desc2: y,x,dim (height2 by width2 by desc_dim)
// mask1: y,x (height1 by width1, 0:bg, 1:fg, 2:NA
// mask2: y,x (height2 by width2)
// dist_mat: y,x,dy,dx (grid range, but only selected disparities)
{
	//double start = omp_get_wtime();
	
	float* desc1 = descriptor1.desc;
	int width1 = descriptor1.width;
	int height1 = descriptor1.height;
	int desc_dim = descriptor1.desc_dim;

	float* desc2 = descriptor2.desc;
	int width2 = descriptor2.width;
	int height2 = descriptor2.height;
	
	ImSize im_size1, im_size2;
	im_size1.height = height1;
	im_size1.width = width1;
	im_size2.height = height2;
	im_size2.width = width2;

	int n_pixel1 = width1*height1;
	int n_pixel2 = width2*height2;
	int search_size_x = params.search_size_x;
	int search_size_y = params.search_size_y;

	int dist_mat_width(0), dist_mat_height(0);
	int disp_x_range[2], disp_y_range[2];
	int dist_mat_dx_range[2], dist_mat_dy_range[2];
	int n_disp_x(0), n_disp_y(0);

	ComputeDistMatSize(im_size1, im_size2,
		               grid_size,
					   search_size_x, search_size_y,
					   dist_mat_width, dist_mat_height,
					   dist_mat_dx_range, dist_mat_dy_range,
					   disp_x_range, disp_y_range,
					   n_disp_x, n_disp_y);

	int n_dist_mat_xy_elem = dist_mat_width*dist_mat_height;
	int n_dist_mat_dx = (dist_mat_dx_range[1] - dist_mat_dx_range[0])/grid_size + 1;
	int n_dist_mat_dy = (dist_mat_dy_range[1] - dist_mat_dy_range[0])/grid_size + 1;
	int n_disp = n_dist_mat_dx * n_dist_mat_dy;
			
#if defined __OPEN_MP
	#pragma omp parallel for	
#endif

	for ( int pos = 0; pos < n_dist_mat_xy_elem; pos++) {
		int x_idx = pos/dist_mat_height;
		int y_idx = pos%dist_mat_height;
		int x1 = x_idx*grid_size;
		int y1 = y_idx*grid_size;
		int ptr1_desc_idx = y1 + x1*height1;

		// for each disparity
		for (int i=0 ; i < n_disp ; i++) {
			
			int dx = disparity[2*i];
			int dy = disparity[2*i+1];

			int x2 = x1 + dx;
			if ( x2 < 0 || x2 >= width2 ) continue;

			int y2 = y1 + dy;
			if ( y2 < 0 || y2 >= height2 ) continue;

			// dist_mat idx at (y1,x1,dy,dx) -- grid range
			int dist_mat_idx = pos + i*n_dist_mat_xy_elem;
			
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