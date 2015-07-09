#pragma once

#include <vector>
#include <bitset>
#include <utility>

using namespace std;

struct MatchParams {
	int sample_grid_size;
	int n_tree_level;
	int search_size_x;
	int search_size_y;
	float truncate_const;	
	float mask_mismatch_const;
};

struct Descriptor{
	float* desc;
	int width;
	int height;
	int desc_dim;
};

struct BinaryDescriptor{
	int width;
	int height;
	static const int desc_dim = 256;
	bitset<desc_dim>* desc;
};

struct ImSize{
	int width;
	int height;
};

struct Disparity{
	int dx;
	int dy;
	float overlap;
	float score;
};

struct Point2D {
	int x;
	int y;
};

class SPNodeMatch {

public:
    
	void NodeMatch(Descriptor descriptor1, 
				   Descriptor descriptor2,
			       MatchParams& params, 
				   float* node_cost,
				   int* node_disparity);

	void NodeMatch(Descriptor descriptor1, 
				   Descriptor descriptor2,
				   unsigned char* mask1,
				   unsigned char* mask2,
				   MatchParams& params,
				   float* node_cost,
				   int* node_disparity);

	void NodeMatch(Descriptor descriptor1, 
				   Descriptor descriptor2,
			       MatchParams& params, 
				   float scale,
			       float* node_cost,
				   int* node_disparity);

	void ComputeDistMatSize(ImSize im_size1, 
							ImSize im_size2,
				     		int grid_size, 
							int search_size_x, int search_size_y,
							int& dist_mat_width, int& dist_mat_height,
							int dist_mat_dx_range[2], int dist_mat_dy_range[2],
							int disp_x_range[2], int disp_y_range[2],
							int &n_disp_x, int &n_disp_y);

	void ComputeDescDist(Descriptor descriptor1, Descriptor descriptor2, 
						 int grid_size,
						 MatchParams& params,
						 float* dist_mat,
						 int* disparity);

	void ComputeDescDist(Descriptor descriptor1, Descriptor descriptor2, 
						 unsigned char* mask1, unsigned char* mask2,
						 int grid_size,
						 MatchParams& params,
						 float* dist_mat,
						 int* disparity);

	void ComputeDescDist(Descriptor descriptor1, Descriptor descriptor2, 
						 int grid_size, float scale,
						 MatchParams& params,
						 float* dist_mat,
						 int* disparity);

	void ComputeNodeCost(float* dist_mat, 
					     ImSize im_size1, ImSize im_size2,
					     MatchParams& params, 
					     int grid_size,
					     float* node_cost);
		                  
					   
};
