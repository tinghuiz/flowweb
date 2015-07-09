#include "SPBP.h"
#include "DT.h"
#include "mex.h"

#include <cmath>
#include <cfloat>
#include <string.h>

#include "omp.h"
#include "openmpflag.h"

SPBP::SPBP(SPBPParams& params)
: bp_params(params), sp_graph(params.n_sp_level)
{
	int n_state = bp_params.n_state;
	int n_node = sp_graph.GetNumNode();
	int n_edge = sp_graph.GetNumEdge();

	// initialize by dummy value
	m_deformation_coeff = 0.0f;
	m_truncate_value = 0.0f;

	// allocate memory
	message = new float[n_edge*n_state*2]; // 2--> incoming/outgoing
	belief = new float[n_node*n_state];
	unary = new float[n_node*n_state];
	optimal_state = new int[n_node];
}

SPBP::~SPBP()
{
	delete [] message;
	delete [] belief;
	delete [] unary;
	delete [] optimal_state;
}


void SPBP::RunBP(float* unary_potential, float deformation_coeff, float truncate_value, int max_iter)
{
	// parameters
	m_deformation_coeff = deformation_coeff;
	m_truncate_value = truncate_value;

	// initialize belief and message
	int n_state = bp_params.n_state;
	int n_node = sp_graph.GetNumNode();
	int n_edge = sp_graph.GetNumEdge();
	
	memcpy(unary, unary_potential, sizeof(float)*n_state*n_node);
	memcpy(belief, unary_potential, sizeof(float)*n_state*n_node);
	for ( int i =0; i < n_edge*n_state*2 ; i++) {
		message[i] = 0.0f;
	}

	// message passing
	for ( int i=0 ; i < max_iter ; i++) {
		TopdownMessagePassing();
		BottomupMessagePassing();
		ForwardMessagePassing();
		BackwardMessagePassing();
		bool change = ComputeOptimalState();
		if ( !change ) {
			mexPrintf("iter: %d\n", i+1);
			break;
		}
	}
}

void SPBP::UpdateMessage(int src_idx, int dst_idx)
{

	if ( src_idx == -1 || dst_idx == -1 ) {
		return;
	}
	
	// parameters
	int n_level = bp_params.n_sp_level;
	int n_state = bp_params.n_state;
	int n_xstate = bp_params.n_xstate;
	int n_ystate = bp_params.n_ystate;

	// message from dst to src
	int edge_dst2src_idx = sp_graph.GetEdgeIdx(dst_idx, src_idx);
	float* dst2src_message = message + edge_dst2src_idx*n_state;

	// belief of dst
	int dst_node_offset = dst_idx*n_state;
	float* dst_belief = belief + dst_node_offset;
	// subtract a message from src --- for update belief 
	int edge_src2dst_idx = sp_graph.GetEdgeIdx(src_idx, dst_idx);
	float* src2dst_message = message + edge_src2dst_idx*n_state;
	for ( int kk=0 ; kk < n_state; kk++) { 
		dst_belief[kk] -= src2dst_message[kk];
	}

	// belief of src
	int src_node_offset = src_idx*n_state;
	float* src_belief = belief + src_node_offset;

	// subtract a message from dst
	for ( int kk=0 ; kk < n_state; kk++) { 
		src2dst_message[kk] = src_belief[kk] - dst2src_message[kk];
	}
	
	// update the message from src to dst
	DT2d(src2dst_message, n_xstate, n_ystate, m_deformation_coeff, m_truncate_value); 

	// update belief of dst with new message from src (w/ normalization)
	float val(0.0f);
	for ( int kk=0 ; kk < n_state; kk++) { 
		dst_belief[kk] += src2dst_message[kk];
		val += dst_belief[kk];
	}
	val /= n_state;
	// normalize
	for ( int kk=0 ; kk < n_state; kk++) { 
		dst_belief[kk] -= val;
	}	
}

void SPBP::TopdownMessagePassing()
{
	// parameters
	int n_level = bp_params.n_sp_level;
	
	// message passing
	for ( int i=0 ; i < n_level-1; i++) {
		int node_inds_at_i[2];
		sp_graph.GetNodeIndsAtLevelK(i, node_inds_at_i);
		int s_idx = node_inds_at_i[0];
		int e_idx = node_inds_at_i[1];

		// for each node (parent) at level i
		for ( int j= s_idx ; j <= e_idx ; j++) { 
			int p_idx = j;
			int children_inds[4];
			sp_graph.GetChildrenNodeInds(p_idx, children_inds);
			
			// for each child nodes
			for ( int k=0; k < 4 ; k++) { 
				int c_idx = children_inds[k];

				UpdateMessage(p_idx, c_idx);
			}
		}
	}
}

void SPBP::BottomupMessagePassing()
{
	// parameters
	int n_level = bp_params.n_sp_level;
	
	// message passing
	for ( int i=n_level-1 ; i >= 1; i--) {
		
		int node_inds_at_i[2];
		sp_graph.GetNodeIndsAtLevelK(i, node_inds_at_i);
		int s_idx = node_inds_at_i[0];
		int e_idx = node_inds_at_i[1];

		// for each node at level i
		for ( int j= s_idx ; j <= e_idx ; j++) { 
			int c_idx = j;
			int p_idx = sp_graph.GetParentNodeIdx(c_idx);
			
			UpdateMessage(c_idx, p_idx);

		}
	}
}

void SPBP::ForwardMessagePassing()
{
	// parameters
	int n_level = bp_params.n_sp_level;
	
	for ( int i= 1 ; i < n_level ; i++) {

		int node_inds_at_i[2];
		sp_graph.GetNodeIndsAtLevelK(i, node_inds_at_i);
		int offset_idx = node_inds_at_i[0]; 
		int n_bin = (int)pow(2.0,i);

		int max_x = 2*(n_bin-1);
		for ( int x = 0; x <= max_x ; x++) {
			
			for ( int x1 = x; x1 >= 0 ; x1-- ) {
				int y1 = x - x1;

				if ( x1 < 0 || x1 >= n_bin || y1 < 0 || y1 >= n_bin ) 
					continue;

				int src_idx = offset_idx + y1 + x1*n_bin;
				// to the right
				if ( x1+1 < n_bin ) {
					int right_idx = offset_idx + y1 + (x1+1)*n_bin;
					UpdateMessage(src_idx, right_idx);
				}
				// to the down
				if  (y1+1 < n_bin ) {
					int down_idx = offset_idx + y1+1 + x1*n_bin;
					UpdateMessage(src_idx, down_idx);
				}
			}
		}
	}
}

void SPBP::BackwardMessagePassing()
{
	// parameters
	int n_level = bp_params.n_sp_level;
	
	for ( int i= 1 ; i < n_level ; i++) {

		int node_inds_at_i[2];
		sp_graph.GetNodeIndsAtLevelK(i, node_inds_at_i);
		int offset_idx = node_inds_at_i[0]; // offset
		int n_bin = (int)pow(2.0,i);

		int max_x = 2*(n_bin-1);
		for ( int x = max_x ; x >=0 ; x--) {
			
			for ( int x1 = x; x1 >= 0 ; x1-- ) {
				int y1 = x - x1;

				if ( x1 < 0 || x1 >= n_bin || y1 < 0 || y1 >= n_bin ) 
					continue;

				int src_idx = offset_idx + y1 + x1*n_bin;
				
				// to the left
				if ( x1-1 >= 0 ) {
					int left_idx = offset_idx + y1 + (x1-1)*n_bin;
					UpdateMessage(src_idx, left_idx);
				}
				// to the up
				if  (y1-1 >= 0 ) {
					int up_idx = offset_idx + y1-1 + x1*n_bin;
					UpdateMessage(src_idx, up_idx);
				}
			}
		}
	}
}


bool SPBP::ComputeOptimalState()
{
	int n_node = sp_graph.GetNumNode();
	int n_state = bp_params.n_state;

	int* optimal_state_ = new int[n_node];

	for (int i=0 ; i < n_node ; i++) {
		int node_offset = i*n_state;
		float* belief_i = belief + node_offset;
	
		float min_val(FLT_MAX);
		int min_idx(-1);	
		for (int j=0 ; j < n_state; j++) {
			if ( belief_i[j] < min_val ) {
				min_val = belief_i[j];
				min_idx = j;
			}
		}
		optimal_state_[i] = min_idx;
	}

	// check if states are updated
	int cnt(0);
	for ( int i=0 ; i < n_node ; i++) {
		if ( optimal_state[i] != optimal_state_[i] ) break;

		cnt++;
	}
	bool change(true);
	if  ( cnt == n_node ) change = false;

	memcpy(optimal_state, optimal_state_, sizeof(int)*n_node);

	return change;
}

void SPBP::GetOptimalState(int* optimal_state_)
{
	int n_node = sp_graph.GetNumNode();
	memcpy(optimal_state_, optimal_state, sizeof(int)*n_node);
}

double SPBP::ComputeEnergy()
{
	// parameter	
	int n_node = sp_graph.GetNumNode();
	int n_edge = sp_graph.GetNumEdge();
	int n_state = bp_params.n_state;
	int n_ystate = bp_params.n_ystate;

	// unary
	double unary_energy(0.0);
	for ( int i=0 ; i < n_node ; i++) {
		int opt_state = optimal_state[i];
		float* unary_i = unary + n_state*i;
		unary_energy += unary_i[opt_state];
	}

	// pairwise
	double pair_energy(0.0);
	for ( int i=0 ; i < n_edge ; i++) {
		pair<int, int> node_inds = sp_graph.GetNodeIdx(i);
		int node1 = node_inds.first;
		int node2 = node_inds.second;
		int state1 = optimal_state[node1];
		int state2 = optimal_state[node2];

		int x1 = state1/n_ystate;
		int y1 = state1%n_ystate;
		int x2 = state2/n_ystate;
		int y2 = state2%n_ystate;

		double dist = abs(double(x1-x2)) + abs(double(y1-y2));
		dist = m_deformation_coeff*min(dist, (double)m_truncate_value);

		pair_energy += dist;
	}

	double energy = unary_energy + pair_energy;
	return energy;

}







