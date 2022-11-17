#include "unispg_A.h"
// #include "GBitVec.h"
// #include <float.h>
// #include <limits.h>

// #include <iostream>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"

void UnispgGp_APPLY::AddUnispg(int s, UnispgGp_APPLY* unispg) {
	int new_refstart = 0;
	int new_refend = 0;
	if (this->graph_num[0] == 0 && this->graph_num[1] == 0) {
		new_refstart = unispg->refstart;// the start of the first node.
		new_refend = unispg->refend;; // the end of the last node.
	} else {
		new_refstart = this->refstart < unispg->refstart ?  this->refstart : unispg->refstart;// the start of the first node.
		new_refend = this->refend > unispg->refend ?  this->refend : unispg->refend;; // the end of the last node.
	}

	// fprintf(stderr, "## new_refstart: %d; new_refend: %d\n", new_refstart, new_refend);

	this->refstart = new_refstart;
	this->refend = new_refend;
	int g_idx = this->graph_num[s];

	// fprintf(stderr, ">> samples.Count(): %d\n", this->samples.Count());
	// fprintf(stderr, "**** 'Copy unispg->no2gnode_unispg to the unispg' check !!!!\n");

	// fprintf(stderr, ">> 'Before copy: 'unispg->no2gnode_unispg[%d][%d]: %d\n", s, 0, unispg->no2gnode_unispg[s][0].Count());	
	int copy_target_g_num = unispg->no2gnode_unispg[s][0].Count();
	// fprintf(stderr, ">> 'Before copy: 'this->no2gnode_unispg[%d][%d]: %d\n", s, g_idx, this->no2gnode_unispg[s][g_idx].Count());
	if (copy_target_g_num > 0) {
		this->no2gnode_unispg[s][g_idx] = new GPVec<CGraphnodeUnispg>(unispg->no2gnode_unispg[s][0]);
		// fprintf(stderr, ">> 'After copy: 'unispg->no2gnode_unispg[%d][%d]: %d\n", s, 0, unispg->no2gnode_unispg[s][0].Count());
		// fprintf(stderr, ">> 'After copy: 'this->no2gnode_unispg[%d][%d]: %d\n", s, g_idx, this->no2gnode_unispg[s][g_idx].Count());
	}
	/*{ // DEBUG
		fprintf(stderr, "**** Post 'Copy_new_no2gnode_unispg_2_no2gnode_unispg' check !!!!\n");
		for(int s=0;s<2;s++) {
			fprintf(stderr, ">> 'Before copy: 'unispg->no2gnode_unispg[%d][%d]: %d\n", s, 0, unispg->no2gnode_unispg[s][0].Count());	
			int copy_target_g_num = unispg->no2gnode_unispg[s][0].Count();
			int g_idx = this->no2gnode_unispg[s]->Count();
			if (copy_target_g_num > 0) {
				for (int gg=0; gg<g_idx; gg++) {
					for (int n=0; n<this->no2gnode_unispg[s][gg].Count(); n++) {
						fprintf(stderr, "\t>> this->no2gnode_unispg[%d][%d][%d]: %d\n", s, gg, n, this->no2gnode_unispg[s][gg].Get(n)->nodeid);
					}
				}
			}
		};
	}
	*/

	// adding unispg->node_nums[s][i] into this->node_nums[s] // how many nodes are in a certain graph g, on strand s: graphno[s][g]
	int node_num = unispg->node_nums[s][0];
	this->node_nums[s].cAdd(node_num);

	// adding unispg->edge_nums[s][i] into this->edge_nums[s] // how many edges are in a certain graph g, on strand s: edge_nums[s][g]
	int edge_num = unispg->edge_nums[s][0];
	this->edge_nums[s].cAdd(edge_num);


	// fprintf(stderr, "unispg->node_nums[s][0]: %d\n", unispg->node_nums[s][0]);
	// fprintf(stderr, "unispg->edge_nums[s][0]: %d\n", unispg->edge_nums[s][0]);
	// fprintf(stderr, "unispg->gpos[s][0].Count(): %d\n", unispg->gpos[s][0].Count());
	
	unispg->gpos[s][0].startIterate();
	int* key;
	int idx=node_num;
	int edge_num_lcl = 1;
	while (true) {
		// fprintf(stderr, "value: %d;\n", val);
		key = unispg->gpos[s][0].Next(idx);
		if (key == NULL) break;
		// fprintf(stderr, ">> counter: %d\n", counter);
		// fprintf(stderr, ">> idx: %d  (g_idx: %d)\n", idx, g_idx);
		// fprintf(stderr, ">> key: %d; idx: %d;\n", *key, idx);

		
		int* pos = this->gpos[s][g_idx][*key];
		if(pos==NULL) {
			this->gpos[s][g_idx].Add(*key, idx);
			pos = this->gpos[s][g_idx][*key];
			// fprintf(stderr, "\t>>> (1) key: %d; *pos: %d\n", *key, *pos);
		} else {
			// fprintf(stderr, "\t>>> (2) key: %d; *pos: %d\n", *key, *pos);
		}
		edge_num_lcl += 1;
	}

	// this->gpos[s][g_idx].startIterate();
	// int* key_test;
	// int idx_test=unispg->node_nums[s][0];
	// while (true) {
	// 	key_test = this->gpos[s][g_idx].Next(idx_test);
	// 	if (key_test == NULL) break;
	// 	// fprintf(stderr, ">> counter: %d\n", counter);
	// 	fprintf(stderr, ">> idx_test: %d\n", idx_test);
	// 	fprintf(stderr, ">> key_test: %d; idx_test: %d;\n", *key_test, idx_test);
	// 	// counter += 1;
	// 	int *val = this->gpos[s][g_idx][*key_test];
	// 	fprintf(stderr, ">> val: %d\n", *val);
	// }

	this->graph_num[s] += 1;
}

void UnispgGp_APPLY::set_N_E_num(int s, int node_num, int edge_num) {
	// fprintf(stderr, "Adding node_num: %d;  edge_num: %d \n", node_num, edge_num);
    this->node_nums[s].cAdd(node_num);
    this->edge_nums[s].cAdd(edge_num);
}

#endif