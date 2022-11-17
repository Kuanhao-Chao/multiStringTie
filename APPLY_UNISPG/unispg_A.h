#ifndef __UNISPG_A_H__
#define __UNISPG_A_H__
#pragma once

#include "global_params_A.h"
#include "../unispg.h"

#include <unordered_map>

struct UnispgGp_APPLY:public UnispgGp {
	public:
		int s_single_dot = -1;
    	GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    	GVec<int> lastgpos[2];

		// GVec<int> node_nums[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
		// GVec<int> edge_nums[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]

		UnispgGp_APPLY(int refstart_i, int refend_i, int g_num, GStr refseq_i) {
			refstart = refstart_i;
			refend = refend_i;
			refseq = refseq_i;
			for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
				int s=sno/2; // adjusted strasnd due to ignoring neutral strand
				no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[g_num];
				gpos[s] = new GIntHash<int>[g_num];
			}
		}

		UnispgGp_APPLY() {
			for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
				int s=sno/2; // adjusted strand due to ignoring neutral strand
				no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
				gpos[s] = new GIntHash<int>[20000];			
			}
		}

		void Clear() {
			// fprintf(stderr, "**** Start Clearing !!!! \n ");
			refstart = 0;
			refend = 0;
			refseq = "";
			for(int i=0;i<2;i++) {
				if (graph_num[i]) {
					// delete [] no2gnode_unispg[i];
					no2gnode_unispg[i]->Clear();
					no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[200];
					node_nums[i].Clear();
					edge_nums[i].Clear();
					graph_num[i] = 0;
					// delete [] gpos[i];
					gpos[i]->Clear();
					gpos[i] = new GIntHash<int>[200];	
				}
			};
		}

		void set_N_E_num(int s, int node_num, int edge_num);
		void AddUnispg(int s, UnispgGp_APPLY* unispg);
};
#endif

