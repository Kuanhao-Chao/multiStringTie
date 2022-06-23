#include "GArgs.h"
#include "rlink.h"

bool segs_overlap(int s1_start, int s1_end, int s2_start, int s2_end);
void calculate_ovp_coverage(int pos_start, int pos_end, int neg_start, int neg_end, int refstart, int refend, float bundle_coverage_ratio_pos, float bundle_coverage_ratio_neg, float& pos_cov, float& neg_cov, GVec<float>* bpcov);