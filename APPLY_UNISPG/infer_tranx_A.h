#ifndef __INFER_TRANX_A_H__
#define __INFER_TRANX_A_H__

#include "global_params_A.h"
#include "tree_pattern_A.h"
#include "parse_reads_A.h"
#include "processTransfrags_A.h"
#include "findTranscripts_A.h"
#include "rlink_A.h"

void infer_transcripts_APPLY_UNISPG(BundleData* bundle, UnispgGp_APPLY* unispgs);

void create_graph_param(int s, int g, UnispgGp_APPLY* unispgs, GPVec<CTransfrag> transfrag, GIntHash<int>** gpos, int& lastgpos);

void graph_dfs(int s, int g, UnispgGp_APPLY* unispgs, CGraphnodeUnispg* node, GVec<bool>& visit);

GBitVec traverse_dfs_APPLY_UNISPG(int s,int g,CGraphnodeUnispg *node,CGraphnodeUnispg *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnodeUnispg> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos);
#endif