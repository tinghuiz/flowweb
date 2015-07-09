#include "SPGraph.h"
#include <cmath>
#include <string.h>

SPGraph::SPGraph(int n_level)
: n_level(n_level), node2level(NULL), node2edge(NULL)
{
	n_node = ((int)pow(4.0,n_level)-1)/3;
	// node2level 
	node2level = new int[n_node];
	n_edge = n_node - 1; // initialized by inter-level edges
	
	for (int i =0; i < n_level ; i++) {
		int offset = ((int)pow(4.0,i)-1)/3;
		int n_node_level_i = (int)pow(4.0,i);
		int n_bin_level_i = (int)pow(2.0,i);

		for ( int j=0 ; j < n_node_level_i ; j++) {
			int node_idx = offset + j;
			node2level[node_idx] = i;
		}
		n_edge += 2*n_bin_level_i*(n_bin_level_i-1); // add by intra-level edges
	}
	edge2node.resize(n_edge);


	// node2edge, edge2node
	// inter-level edge
	node2edge = new int[n_node*n_node];
	memset(node2edge, -1, sizeof(int)*n_node*n_node);
	int edge_idx(0);
	for (int i=0 ; i < n_level; i++) {
		int offset = ((int)pow(4.0,i)-1)/3;
		int n_node_level_i = (int)pow(4.0,i);

		for ( int j=0 ; j < n_node_level_i ; j++) {
			int children_inds[4];
			int node_idx = offset + j;
			GetChildrenNodeInds(node_idx, children_inds);
			for ( int k=0 ; k < 4 ; k++) {
				int c_idx = children_inds[k];
				if ( c_idx >= 0 ) {
					int idx1 = node_idx + c_idx*n_node;
					int idx2 = c_idx + node_idx*n_node;
					if ( node2edge[idx1] == -1 && node2edge[idx2] == -1 ) {
						node2edge[idx1] = 2*edge_idx;
						node2edge[idx2] = 2*edge_idx+1;
						edge2node[edge_idx].first = idx1;
						edge2node[edge_idx].second = idx2;
						edge_idx++;
					}
				}
			}
		}
	}
	// intra-level edge
	for (int i=0 ; i < n_level; i++) {
		int offset = ((int)pow(4.0,i)-1)/3;
		int n_node_level_i = (int)pow(4.0,i);

		for ( int j=0 ; j < n_node_level_i ; j++) {
			int neighbor_inds[4];
			int node_idx = offset + j;
			GetNeighborNodeInds(node_idx, neighbor_inds);
			for ( int k=0 ; k < 4 ; k++) {
				int n_idx = neighbor_inds[k];
				if ( n_idx >= 0 ) {
					int idx1 = node_idx + n_idx*n_node;
					int idx2 = n_idx + node_idx*n_node;
					if ( node2edge[idx1] == -1 && node2edge[idx2] == -1 ) {
						node2edge[idx1] = 2*edge_idx;
						node2edge[idx2] = 2*edge_idx+1;
						edge2node[edge_idx].first = idx1;
						edge2node[edge_idx].second = idx2;
						edge_idx++;
					}
				}
			}
		}
	}
	//assert(n_edge == edge_idx);
}

SPGraph::~SPGraph()
{
	delete [] node2level;
	delete [] node2edge;
}


SPGraph::SPGraph(const SPGraph& that)
: n_level(that.n_level), n_node(that.n_node), n_edge(that.n_edge)
{
	node2level = new int[n_node];
	memcpy(node2level, that.node2level, sizeof(int)*n_node);
}


SPGraph& SPGraph::operator =(const SPGraph& that)
{
	if ( this == &that) return *this;
	delete [] node2level;

	n_node = that.n_node;
	n_edge = that.n_edge;
	n_level = that.n_level;

	node2level = new int[n_node];
	memcpy(node2level, that.node2level, sizeof(int)*n_node);

	return *this;
}

int SPGraph::GetParentNodeIdx(int node_idx)
{
	int level_idx = node2level[node_idx];
	
	int p_idx = -1;
	if ( level_idx > 0 ) {
		int offset = ((int)pow(4.0,level_idx)-1)/3;
		int n_bin = (int)pow(2.0,level_idx);
		int x_idx = (node_idx - offset)/n_bin;
		int y_idx = (node_idx - offset)%n_bin;

		int p_x_idx = x_idx/2;
		int p_y_idx = y_idx/2;
		int p_n_bin = n_bin/2;
		int p_offset = ((int)pow(4.0,level_idx-1)-1)/3;
		int p_idx = p_offset + p_y_idx + p_x_idx*p_n_bin;
	}

	return p_idx;
}
void SPGraph::GetChildrenNodeInds(int node_idx, int children_inds[4])
{
	int level_idx = node2level[node_idx];

	memset(children_inds, -1, sizeof(int)*4);

	if ( level_idx < n_level-1 ) {
		int offset = ((int)pow(4.0,level_idx)-1)/3;
		int n_bin = (int)pow(2.0,level_idx);
		int x_idx = (node_idx - offset)/n_bin;
		int y_idx = (node_idx - offset)%n_bin;

		int c_offset = ((int)pow(4.0,level_idx+1)-1)/3;

		for (int x=0; x <= 1; x++) {
			for (int y=0; y <= 1; y++) {
				int c_x_idx = 2*x_idx + x;
				int c_y_idx = 2*y_idx + y;
				int c_n_bin = 2*n_bin;
				
				int c_idx = c_offset + c_y_idx + c_x_idx*c_n_bin;
				children_inds[y+x*2] = c_idx;
			}
		}
	}
}

void SPGraph::GetNeighborNodeInds(int node_idx, int neighbor_inds[4])
{
	int level_idx = node2level[node_idx];

	memset(neighbor_inds, -1, sizeof(int)*4);

	int offset = ((int)pow(4.0,level_idx)-1)/3;
	int n_bin = (int)pow(2.0,level_idx);
	int x_idx = (node_idx - offset)/n_bin;
	int y_idx = (node_idx - offset)%n_bin;

	int cnt(0);

	int x[4];
	int y[4];
	x[0] = -1; y[0] = 0; // left
	x[1] = 1; y[1] = 0; // right
	x[2] = 0; y[2] = -1; // up
	x[3] = 0; y[3] = 1; // down
	
	for ( int i=0; i < 4 ; i++ ) {
	
		int n_x_idx = x_idx + x[i];
		int n_y_idx = y_idx + y[i];
		if ( n_x_idx >=0 && n_y_idx >=0 && 
			 n_x_idx < n_bin && n_y_idx < n_bin) {
			int n_idx = offset + n_y_idx + n_x_idx*n_bin;
			neighbor_inds[i] = n_idx;
		}
	}

}

void SPGraph::GetNodeIndsAtLevelK(int level_idx, int node_inds[2])
{
	int offset = ((int)pow(4.0,level_idx)-1)/3;
	int n_node_level_i = (int)pow(4.0,level_idx);
	node_inds[0] = offset;
	node_inds[1] = offset + n_node_level_i - 1;
}

int SPGraph::GetEdgeIdx(int node_idx1, int node_idx2)
{
	int idx = node_idx2 + node_idx1*n_node;
	return node2edge[idx];
}