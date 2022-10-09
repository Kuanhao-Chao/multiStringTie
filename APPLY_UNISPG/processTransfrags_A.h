#ifndef __PROCESSTRANSFRAGS_A_H__
#define __PROCESSTRANSFRAGS_A_H__

#include "global_params_A.h"
#include "../unispg.h"
#include "../helper.h"

void process_transfrags_APPLY_UNISPG(int s, int gno,int edgeno,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no, GIntHash<int> &gpos, GList<CPrediction>& pred);

int compatible_long_APPLY_UNISPG(int* t,int *len,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,int gno,GIntHash<int> &gpos);

bool assign_incomplete_trf_to_nodes_APPLY_UNISPG(int t,int n1, int n2,GPVec<CGraphnodeUnispg>& no2gnode);

bool trf_real_APPLY_UNISPG(int t,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag, GIntHash<int> &gpos,int gno);

// returns true if at list one node gets assigned the incomplete transfrag
bool binary_insert_trf_to_node_APPLY_UNISPG(int t, GVec<int>& trf,int first,int last);
#endif