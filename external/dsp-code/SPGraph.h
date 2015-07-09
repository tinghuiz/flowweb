#pragma once

#include <vector>
#include <utility>

using namespace std;

class SPGraph {
public:
	explicit SPGraph(int n_level);
	SPGraph(const SPGraph& that);
	SPGraph& operator = ( const SPGraph& that);

	virtual ~SPGraph();

public:
	void GetNodeIndsAtLevelK(int level_idx, int node_inds[2]);
	int GetParentNodeIdx(int node_idx);
	void GetChildrenNodeInds(int node_idx, int children_inds[4]);
	void GetNeighborNodeInds(int node_idx, int neighbor_inds[4]);
	int GetNumEdge()const { return n_edge;};
	int GetNumNode()const { return n_node;};
	int GetEdgeIdx(int node_idx1, int node_idx2);
	pair<int, int>& GetNodeIdx(int edge_idx) { return edge2node[edge_idx]; };

private:
	int* node2level;
	int* node2edge;
	vector<pair<int, int> > edge2node;
	int n_node;
	int n_edge;
	int n_level;
};

