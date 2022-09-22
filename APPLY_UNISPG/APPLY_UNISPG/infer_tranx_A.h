#ifndef __INFER_TRANX_A_H__
#define __INFER_TRANX_A_H__

#include "global_params_A.h"
#include "tree_pattern_A.h"
#include "parse_reads_A.h"
#include "processTransfrags_A.h"
#include "rlink_A.h"

void infer_transcripts_APPLY_UNISPG(BundleData* bundle, GPVec<UnispgGp_APPLY>** graphs_vec);

void create_graph_param(int s, int g, GPVec<UnispgGp_APPLY>** graphs_vec, GPVec<CTransfrag> transfrag, GIntHash<int>** gpos, int& lastgpos);

void graph_dfs(int s, int g, GPVec<UnispgGp_APPLY>** graphs_vec, CGraphnodeUnispg* node, GVec<bool>& visit);

#endif