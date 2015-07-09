#pragma once

#include <vector>
#include <utility>

//#include "SPNodeMatch.h"

struct PixelMatchParams {
	int search_grid_size;
	int search_size_x;
	int search_size_y;
	float deformation_coeff;
	float truncate_const;
	float mask_mismatch_const;
};

struct Descriptor{
	float* desc;
	float* norm;
	int width;
	int height;
	int desc_dim;
};


class PixelMatch{

public:
	
	void RunPixelMatch(Descriptor descriptor1, 
					   Descriptor descriptor2,
					   int* pixel_disparity, 
					   float* match_cost,
				       PixelMatchParams& params);

	void RunPixelMatch(Descriptor descriptor1, 
					   Descriptor descriptor2,
					   int* bbox,
					   int* bbox_disparity,
					   int* pixel_disparity,
					   float* match_cost,
					   int n_bbox,
				       PixelMatchParams& params);

	void RunPixelMatch(Descriptor descriptor1, 
					   Descriptor descriptor2,
					   unsigned char* mask1,
					   unsigned char* mask2,
					   int* pixel_disparity, 
					   float* match_cost,
				       PixelMatchParams& params);


private:

	void ComputeDescDist(Descriptor descriptor1, 
						 Descriptor descriptor2, 
						 PixelMatchParams& params,
						 float* dist_mat,
						 int* pixel_disparity,
						 int* search_range,
						 int n_search);

	void ComputeDescDist(Descriptor descriptor1, 
						 Descriptor descriptor2, 
						 unsigned char* mask1,
						 unsigned char* mask2,
						 PixelMatchParams& params,
						 float* dist_mat,
						 int* pixel_disparity,
						 int* search_range,
						 int n_search);

	void ComputeBestDisparity(float* dist_mat, 
							  int* pixel_disparity,
							  float* match_cost,
							  int* search_range,
							  float deformation_coeff,
					          int width, int height, int n_search);

	void UpdateDescDist(Descriptor descriptor1,
						Descriptor descriptor2,
						PixelMatchParams& params,
						float* match_cost,
						int* pixel_disparity,
						int bbox[4],
						int ref_disparity[2],						
						int* search_range,
						float* deformation_cost,
						int n_search);
	
};

