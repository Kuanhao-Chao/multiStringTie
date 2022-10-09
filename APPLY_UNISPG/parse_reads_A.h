#ifndef __PARSE_READ_A_H__
#define __PARSE_READ_A_H__

#include "global_params_A.h"
#include "../helper.h"
#include "unispg_A.h"
#include "rlink_A.h"

void get_fragment_pattern_APPLY_UNISPG(BundleData* bundle, GList<CReadAln>& readlist, int n, int np, float readcov, GPVec<UnispgGp_APPLY>** graphs_vec, int* global_gidx, GPVec<CTransfrag> **transfrag, GPVec<CMTransfrag>** mgt, CTreePat ***tr2no);

void get_read_pattern_APPLY_UNISPG(int s, float readcov, float rprop, GList<CReadAln>& readlist, int n, GVec<int>& nodes, GPVec<UnispgGp_APPLY>** graphs_vec, GBitVec& pat, int* global_gidx, GPVec<CTransfrag> **transfrag, GPVec<CMTransfrag>** mgt, CTreePat ***tr2no);

CTransfrag* update_read_pattern_abund_APPLY_UNISPG(int s,int g, int node_num, GIntHash<int>& gpos, GBitVec& pat, float abundance, GVec<int>& nodes, GPVec<CTransfrag> **transfrag, GPVec<CMTransfrag>** mgt, CTreePat ***tr2no);
#endif