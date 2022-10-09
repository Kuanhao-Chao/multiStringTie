#ifndef __UNISPG_A_H__
#define __UNISPG_A_H__
#pragma once

#include "global_params_A.h"
#include "../unispg.h"

#include <unordered_map>

struct UnispgGp_APPLY:public UnispgGp {
	public:
		int graph_num[2] = {0};  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    	GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    	GVec<int> lastgpos[2];

		// GVec<int> node_nums[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
		// GVec<int> edge_nums[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]

		UnispgGp_APPLY(int refstart_i, int refend_i, int g_num) {
			refstart = refstart_i;
			refend = refend_i;
			for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
				int s=sno/2; // adjusted strasnd due to ignoring neutral strand
				no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[g_num];
				gpos[s] = new GIntHash<int>[g_num];
			}
		}
		void set_N_E_num(int s, int node_num, int edge_num);
		void AddUnispg(int s, UnispgGp_APPLY* unispg);
};
#endif

