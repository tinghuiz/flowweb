#pragma once

#include <vector>
#include <utility>

#include "SPGraph.h"

using namespace std;

struct SPBPParams {
	int n_sp_level;
	int n_state;
	int n_xstate;
	int n_ystate;
};



class SPBP {
public:
	SPBP(SPBPParams& params);
	virtual ~SPBP();
	
public:

	void RunBP(float* unary_potential, float deformation_coeff, float truncate_value, int max_iter);
	void GetOptimalState(int* optimal_state_);
	double ComputeEnergy();

private:
	
	void TopdownMessagePassing();
	void BottomupMessagePassing();
	void ForwardMessagePassing();
	void BackwardMessagePassing();
	void UpdateMessage(int src_idx, int dst_idx);
	bool ComputeOptimalState();

private:
	float* unary;
	float* belief;
	float* message;
	
	int* optimal_state;
	float m_deformation_coeff;
	float m_truncate_value;
	SPGraph sp_graph;
	SPBPParams bp_params;
};

