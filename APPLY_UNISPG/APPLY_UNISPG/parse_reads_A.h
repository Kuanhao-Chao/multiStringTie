#ifndef __PARSE_READ_A_H__
#define __PARSE_READ_A_H__

#include "global_params_A.h"
#include "../helper.h"
#include "unispg_A.h"

void get_fragment_pattern_APPLY_UNISPG(BundleData* bundle, GList<CReadAln>& readlist, int n, int np, float readcov, GPVec<UnispgGp_APPLY>** graphs_vec, int* global_gidx);

void get_read_pattern_APPLY_UNISPG(int s, float readcov, float rprop, GList<CReadAln>& readlist, int n,GVec<int>& rgno, GVec<int> *rnode, GPVec<UnispgGp_APPLY>** graphs_vec, int* global_gidx);

#endif