#ifndef __HELPER_H__
#define __HELPER_H__

#include "global_params.h"
#include "rlink.h"

bool segs_overlap(int s1_start, int s1_end, int s2_start, int s2_end);

bool segs_overlap_chain(char target_node_strand, int pos_start, int pos_end, int neg_start, int neg_end, int refstart, int refend, float bundle_coverage_ratio_pos, float bundle_coverage_ratio_neg, int& last_ovp_end, char& chaining_hold_strand, float& chaining_hold_cov, float& end_chaining_cov, GVec<float>* bpcov);

void ovp_coverage_push_node(int& g_idx, int& g_num, int& n_idx, int& n_num, bool& reach_end);

void calculate_ovp_coverage(int pos_start, int pos_end, int neg_start, int neg_end, int refstart, int refend, float bundle_coverage_ratio_pos, float bundle_coverage_ratio_neg, float& pos_cov, float& neg_cov, GVec<float>* bpcov);

void redistribute_unstranded_rcov(float* rprop, GVec<float>* bpcov, int refstart, int refend, int rstart, int rend);

int calOverlapLen(uint rstart, uint rend, uint start, uint end);

char* sprintTime();

void addPtFeature(const char* refname, GPtFeature* pf, GArray<GRefPtData>& refpts);

int loadPtFeatures(FILE* f, GArray<GRefPtData>& refpts);
#endif