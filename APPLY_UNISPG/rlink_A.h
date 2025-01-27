#ifndef __RLINK_A_H__
#define __RLINK_A_H__

#pragma once
#include "global_params_A.h"
#include "../helper.h"

void cov_edge_add_APPLY_UNISPG(GVec<float> *bpcov, int sno, int start, int end, float v);

void add_read_to_cov_APPLY_UNISPG(GList<CReadAln>& rd,int n,GVec<float> *bpcov,int refstart, int refend);

void count_good_junctions_APPLY_UNISPG(BundleData* bdata);

CTransfrag *findtrf_in_treepat_APPLY_UNISPG(int gno,GIntHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no);
#endif