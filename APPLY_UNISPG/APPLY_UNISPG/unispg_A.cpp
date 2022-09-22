#include "unispg_A.h"
// #include "GBitVec.h"
// #include <float.h>
// #include <limits.h>

// #include <iostream>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"

void UnispgGp_APPLY::AddUnispg(int s, UnispgGp_APPLY* unispg) {
		// int refstart = 0; // the start of the first node.
		// int refend = 0; // the end of the last node.
		// GPVec<CGraphnodeUnispg>* no2gnode_unispg[2]; // for each graph g, on a strand s, no2gnode_unispg[s][g][i] gives the node i

		// GVec<int> node_nums[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
		// GVec<int> edge_nums[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]

		// GVec<GStr> samples;
		// CGraphnodeUnispg* source_gp[2];
		// CGraphnodeUnispg* sink_gp[2];

	/*****************************
	 * Update data structure
	 *****************************/   


			// for(int s=0;s<2;s++) {
				
			// 	// fprintf(stderr, ">> new_gidx[%d]: %d\n", s, new_gidx[s]);
			// 	// for (int j=0; j<new_gidx[s]; j++) {
			// 	fprintf(stderr, ">> 'Copy_new_no2gnode_unispg_2_no2gnode_unispg: 'new_no2gnode_unispg[%d][%d]: %d\n", s, 0, new_no2gnode_unispg[s]->Count());
			// 	// GPVec<CGraphnodeUnispg> tmp = new GPVec(new_no2gnode_unispg[s][j]);
			// 	no2gnode_unispg[s][current_gidx[s]] = new GPVec<CGraphnodeUnispg>(new_no2gnode_unispg[s][0]);
			// 	// }
			// 	// no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			// };
			

    // this -> no2gnode_unispg
    // for (int s_itr=0; s_itr<2; s_itr++) {
    //     for (int g_in=0; g_in < unispg->) {

    //     }
    //     this -> no2gnode_unispg[s]

    //     unispg

	// 	node_nums[s];
	// 	edge_nums[s]; 
    // }
    // unispg->no2gnode_unispg

	// /*****************************
	//  * Update the start & end of UnispgGp_APPLY
	//  *****************************/   
    // refstart = refstart_i;
    // refend = refend_i;
}

void UnispgGp_APPLY::set_N_E_num(int s, int node_num, int edge_num) {
	fprintf(stderr, "Adding node_num: %d;  edge_num: %d \n", node_num, edge_num);
    this->node_nums[s].Add(node_num);
    this->edge_nums[s].Add(edge_num);
}

#endif