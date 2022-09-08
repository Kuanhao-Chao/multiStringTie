#include "helper.h"

bool segs_overlap(int s1_start, int s1_end, int s2_start, int s2_end) {
    // fprintf(stderr, "Inside `segs_overlap`\n");
    bool ovp =true;
    // fprintf(stderr, "Process two segments (pre: %d - %d ;  now: %d - %d)!!!\n", s1_start, s1_end, s2_start, s2_end);
    // fprintf(stderr, "\t>>>>> brec: %d - %d\n", s1_start, s1_end);
    if (s1_end <  s2_start) {
        // fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
        ovp = false;
    } else if (s1_start < s2_start && s1_end == s2_start) {
        // fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
        ovp = false;
    } else if (s1_start < s2_start && s1_end > s2_start && s1_end < s2_end) {
        // fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
    } else if (s1_start < s2_start && s1_end > s2_start && s1_end == s2_end) {
        // fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
    } else if (s1_start < s2_start && s1_end > s2_start && s1_end > s2_end) {
        // fprintf(stderr, "\t** Bundle: -----|(s)----------(e)|------\n");
    } else if (s1_start == s2_start && s1_start < s2_end && s1_end < s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
    } else if (s1_start == s2_start && s1_start < s2_end && s1_end == s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
    } else if (s1_start == s2_start && s1_start < s2_end && s1_end > s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end > s2_start && s1_end < s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end == s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end > s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
    } else if (s1_start > s2_start && s1_start == s2_end && s1_end > s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
        ovp = false;
    } else if (s1_start > s2_end) {
        // fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
        ovp = false;
    }
    return ovp;
}


bool segs_overlap_chain(char target_node_strand, int pos_start, int pos_end, int neg_start, int neg_end, int refstart, int refend, float bundle_coverage_ratio_pos, float bundle_coverage_ratio_neg, int& last_ovp_end, char& chaining_hold_strand, float& chaining_hold_cov, float& end_chaining_cov, GVec<float>* bpcov) {
    // fprintf(stderr, "\tInside `segs_overlap`\n");
    // fprintf(stderr, "\t\t Before ~~ (%c) chaining_hold_cov: %f\n", chaining_hold_strand, chaining_hold_cov);
    bool ovp =true;
    // 'chaining_hold_strand' stores the strand of the coverage of the overlapping region.
    // 'chaining_hold_cov' no overlapping of the previous pair => end_strand == '.'


    // float end_chaining_cov = 0;
    float pos_cov_ovp = 0.0;
    float uns_cov_ovp = 0.0;
    float neg_cov_ovp = 0.0;

    fprintf(stderr, "\t>>> last_ovp_end: %d \n", last_ovp_end);


    fprintf(stderr, "\t>> chaining_hold_cov (%c): %f\n", chaining_hold_strand, chaining_hold_cov);

    fprintf(stderr, "\tProcess two segments (pre: %d - %d ;  now: %d - %d)!!!\n", pos_start, pos_end, neg_start, neg_end);
    // fprintf(stderr, "\t>>>>> brec: %d - %d\n", pos_start, pos_end);
    if (pos_end <  neg_start) {
        fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
        ovp = false;
        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            if (chaining_hold_strand == '+' && last_ovp_end > pos_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, pos_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov(%f) + get_cov_sign(2, last_ovp_end[%d]-refstart, pos_end[%d]-refstart, bpcov)(%f) = %f \n", chaining_hold_cov, last_ovp_end, pos_end, get_cov_sign(2, last_ovp_end-refstart, pos_end-refstart, bpcov),end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(2, pos_start[%d]-refstart, pos_end[%d]-refstart, bpcov)(%f) = %f \n", pos_start, pos_end, get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            }
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, neg_end[%d]-refstart, bpcov)(%f) = %f \n", neg_start, neg_end, get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            last_ovp_end = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }
    } else if (pos_start < neg_start && pos_end == neg_start) {
        fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
        ovp = false;
        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            if (chaining_hold_strand == '+' && last_ovp_end > pos_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, pos_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov(%f) +  get_cov_sign(2, last_ovp_end[%d]-refstart, pos_end[%d]-refstart, bpcov)(%f)= %f \n", chaining_hold_cov, last_ovp_end, pos_end, get_cov_sign(2, last_ovp_end-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(2, pos_start[%d]-refstart, pos_end[%d]-refstart, bpcov)(%f) = %f \n", pos_start, pos_end, get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            }
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, neg_end[%d]-refstart, bpcov)(%f) = %f \n", neg_start, neg_end, get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            last_ovp_end = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }
    } else if (pos_start < neg_start && pos_end > neg_start && pos_end < neg_end) {
        fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, neg_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            if (chaining_hold_strand == '+' && last_ovp_end > pos_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov(%f) + get_cov_sign(2, last_ovp_end[%d]-refstart, neg_start[%d]-refstart, bpcov)(%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f \n", chaining_hold_cov, last_ovp_end, neg_start, get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(2, pos_start[%d]-refstart, neg_start[%d]-refstart, bpcov)(%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f)= %f \n", pos_start, neg_start, get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);
            }
            chaining_hold_strand = '-';
            chaining_hold_cov = 0 + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f + neg_cov_ovp(%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", 0, neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg + get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: neg_cov_ovp(%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f](%f) +  get_cov_sign(0, pos_end[%d]-refstart, neg_end[%d]-refstart, bpcov) (%f) = %f \n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, pos_end, neg_end, get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }

    } else if (pos_start < neg_start && pos_end > neg_start && pos_end == neg_end) {
        fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            if (chaining_hold_strand == '+' && last_ovp_end > pos_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov(%f) + get_cov_sign(2, last_ovp_end[%d]-refstart, neg_start[%d]-refstart, bpcov)(%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f \n", chaining_hold_cov, last_ovp_end, neg_start, get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov),  pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(2, pos_start[%d]-refstart, neg_start[%d]-refstart, bpcov)(%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f \n", pos_start, neg_start, get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov),  pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);
            }
            chaining_hold_strand = '-';
            chaining_hold_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            fprintf(stderr, "\t\t** end_chaining_cov: neg_cov_ovp(%f) + uns_cov_ovp*bundle_coverage_ratio_neg(%f) = %f \n", neg_cov_ovp, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);

            // Get the positive `chaining_hold_cov`
            if (chaining_hold_strand == '+') {
                if (last_ovp_end > pos_start) {
                    fprintf(stderr, "\t\t** chaining_hold_cov: chaining_hold_cov (%f) + get_cov_sign(2, last_ovp_end[%d]-refstart, neg_start[%d]-refstart, bpcov) (%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f)\n", chaining_hold_cov, last_ovp_end, neg_start, get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos);
                    chaining_hold_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                } else {
                    fprintf(stderr, "\t\t** chaining_hold_cov: get_cov_sign(2, pos_start[%d]-refstart, neg_start[%d]-refstart, bpcov) (%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f)\n", pos_start, neg_start, get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos);
                    chaining_hold_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                }
            } else {
                chaining_hold_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            }
            chaining_hold_strand = '+';
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }
    } else if (pos_start < neg_start && pos_end > neg_start && pos_end > neg_end) {
        fprintf(stderr, "\t** Bundle: -----|(s)----------(e)|------\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            if (chaining_hold_strand == '+' && last_ovp_end > pos_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos + get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov(%f) + get_cov_sign(2, last_ovp_end[%d]-refstart, neg_start[%d]-refstart, bpcov)(%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) + get_cov_sign(2, neg_end[%d]-refstart, pos_end[%d]-refstart, bpcov)(%f) = %f\n", chaining_hold_cov, last_ovp_end, neg_start, get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, neg_end, pos_end, get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov) + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos + get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov);

                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(2, pos_start[%d]-refstart, neg_start[%d]-refstart, bpcov) (%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) + get_cov_sign(2, neg_end[%d]-refstart, pos_end[%d]-refstart, bpcov) (%f) = %f \n", pos_start, neg_start, get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, neg_end, pos_end, get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            }
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = neg_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            fprintf(stderr, "\t\t** end_chaining_cov: neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);
            if (chaining_hold_strand == '+') {
                if (last_ovp_end > pos_start) {
                    fprintf(stderr, "\t\t** chaining_hold_cov: chaining_hold_cov(%f) + get_cov_sign(2, last_ovp_end[%d]-refstart, neg_start[%d]-refstart, bpcov) (%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f)\n", chaining_hold_cov, last_ovp_end, neg_start, get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos);
                    chaining_hold_cov = chaining_hold_cov + get_cov_sign(2, last_ovp_end-refstart, neg_start-refstart, bpcov)+pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;     
                } else {
                    chaining_hold_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov)+pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                    fprintf(stderr, "\t\t** chaining_hold_cov:  get_cov_sign(2, pos_start[%d]-refstart, neg_start[%d]-refstart, bpcov) (%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", pos_start, neg_start, get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, chaining_hold_cov);
                }
            } else {
                chaining_hold_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov)+pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
                fprintf(stderr, "\t\t** chaining_hold_cov:  get_cov_sign(2, pos_start[%d]-refstart, neg_start[%d]-refstart, bpcov) (%f) + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n",  pos_start, neg_start, get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov), pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, chaining_hold_cov);
            }
            chaining_hold_strand = '+';
            last_ovp_end = neg_end;
        }
    } else if (pos_start == neg_start && pos_start < neg_end && pos_end < neg_end) {
        fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            fprintf(stderr, "\t\t** end_chaining_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);

            chaining_hold_strand = '-';
            chaining_hold_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;

            fprintf(stderr, "\t\t** chaining_hold_cov: neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, chaining_hold_cov);
            last_ovp_end = pos_end;
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg + get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) + get_cov_sign(0, pos_end[%d]-refstart, neg_end[%d]-refstart, bpcov) (%f) = %f\n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, pos_end, neg_end, get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f\n", 0);
            last_ovp_end = pos_end;
        }

    } else if (pos_start == neg_start && pos_start < neg_end && pos_end == neg_end) {
        fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            fprintf(stderr, "\t\t** end_chaining_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);

            chaining_hold_strand = '-';
            chaining_hold_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, chaining_hold_cov);

        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            fprintf(stderr, "\t\t** end_chaining_cov: neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);
            chaining_hold_strand = '+';
            chaining_hold_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, chaining_hold_cov);
        }

    } else if (pos_start == neg_start && pos_start < neg_end && pos_end > neg_end) {
        fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            fprintf(stderr, "\t\t** end_chaining_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) + get_cov_sign(2, neg_end[%d]-refstart, pos_end[%d]-refstart, bpcov) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, neg_end, pos_end, get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov),end_chaining_cov);
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = neg_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            end_chaining_cov = neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            fprintf(stderr, "\t\t** end_chaining_cov: neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);
            chaining_hold_strand = '+';
            chaining_hold_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            last_ovp_end = neg_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, chaining_hold_cov);
        }

    } else if (pos_start > neg_start && pos_start < neg_end && pos_end > neg_start && pos_end < neg_end) {
        fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            fprintf(stderr, "\t\t** end_chaining_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);
            if (chaining_hold_strand == '-') {
                if (last_ovp_end > neg_start) {
                    fprintf(stderr, "\t\t** chaining_hold_cov: chaining_hold_cov (%f) + get_cov_sign(0, last_ovp_end[%d]-refstart, pos_start[%d]-refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f)\n", chaining_hold_cov, last_ovp_end, pos_start, get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg);

                    chaining_hold_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov)+ neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                } else {
                    fprintf(stderr, "\t\t** chaining_hold_cov: %f + get_cov_sign(0, neg_start[%d]-refstart, pos_start[%d]-refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f)\n", 0, neg_start, pos_start, get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg);

                    chaining_hold_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov)+ neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                }
            } else {
                fprintf(stderr, "\t\t** chaining_hold_cov: %f + get_cov_sign(0, neg_start[%d]-refstart, pos_start[%d]-refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f)\n", 0, neg_start, pos_start, get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg);

                chaining_hold_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov)+neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg ;
            }
            chaining_hold_strand = '-';
            last_ovp_end = pos_end;
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            if (chaining_hold_strand == '-' && last_ovp_end > neg_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg + get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);

                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov(%f) + get_cov_sign(0, last_ovp_end[%d]-refstart, pos_start[%d]-refstart, bpcov)(%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) + get_cov_sign(0, pos_end[%d]-refstart, neg_end[%d]-refstart, bpcov)(%f) = %f\n", chaining_hold_cov, last_ovp_end, pos_start, get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, pos_end, neg_end, get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg + get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, pos_start[%d]-refstart, bpcov)(%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) + get_cov_sign(0, pos_end[%d]-refstart, neg_end[%d]-refstart, bpcov)(%f) = %f\n", neg_start, pos_start, get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, pos_end, neg_end, get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            }
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
            last_ovp_end = pos_end;
        }
    } else if (pos_start > neg_start && pos_start < neg_end && pos_end == neg_end) {
        fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            fprintf(stderr, "\t\t** end_chaining_cov: pos_cov_ovp (%f) + uns_cov_ovp*bundle_coverage_ratio_pos (%f) = %f\n", pos_cov_ovp, uns_cov_ovp*bundle_coverage_ratio_pos, end_chaining_cov);
            // Get the negative `chaining_hold_cov`
            if (chaining_hold_strand == '+') {
                chaining_hold_cov =  get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
            } else {
                if (last_ovp_end > neg_start) {
                    fprintf(stderr, "\t\t** chaining_hold_cov: chaining_hold_cov (%f) + get_cov_sign(0, last_ovp_end[%d] -refstart, pos_start[%d] -refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f)\n", chaining_hold_cov, last_ovp_end, pos_start, get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg);

                    chaining_hold_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov)+ neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                } else {
                    fprintf(stderr, "\t\t** chaining_hold_cov: %f + get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp*bundle_coverage_ratio_neg (%f)\n", 0, get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp*bundle_coverage_ratio_neg);

                    chaining_hold_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov)+ neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                }
            }
            chaining_hold_strand = '-';
            last_ovp_end = pos_end;
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            if (chaining_hold_strand == '-' && last_ovp_end > neg_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;

                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov (%f) + get_cov_sign(0, last_ovp_end[%d]-refstart, pos_start[%d]-refstart, bpcov)(%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", chaining_hold_cov, last_ovp_end, pos_start, get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, pos_start[%d]-refstart, bpcov)(%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_start, pos_start, get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);

            }
            chaining_hold_strand = '+';
            chaining_hold_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            last_ovp_end = pos_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }
    } else if (pos_start > neg_start && pos_start < neg_end && pos_end > neg_end) {
        fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
        ovp = true;
        // Calculate the overlapping coverage for pos / uns / neg.
        pos_cov_ovp = get_cov(2, pos_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
        if (pos_cov_ovp <=0.5) pos_cov_ovp = 0;
        if (neg_cov_ovp <=0.5) neg_cov_ovp = 0;
        if (uns_cov_ovp <=0.5) uns_cov_ovp = 0;

        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos + get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) + get_cov_sign(2, neg_end[%d]-refstart, pos_end[%d]-refstart, bpcov) (%f) = %f\n", pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, neg_end, pos_end, get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = neg_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            if (chaining_hold_strand == '-' && last_ovp_end > neg_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov (%f) + get_cov_sign(0, last_ovp_end[%d]-refstart, pos_start[%d]-refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", chaining_hold_cov, last_ovp_end, pos_start, get_cov_sign(0, last_ovp_end-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov) + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, pos_start[%d]-refstart, bpcov) (%f) + neg_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_neg[%f] (%f) = %f\n", neg_start, pos_start, get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov), neg_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_neg, uns_cov_ovp*bundle_coverage_ratio_neg, end_chaining_cov);
            }
            chaining_hold_strand = '+';
            chaining_hold_cov = pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
            last_ovp_end = neg_end;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f + pos_cov_ovp (%f) + uns_cov_ovp[%f]*bundle_coverage_ratio_pos[%f] (%f) = %f\n", 0, pos_cov_ovp, uns_cov_ovp, bundle_coverage_ratio_pos, uns_cov_ovp*bundle_coverage_ratio_pos, chaining_hold_cov);
        }
    } else if (pos_start > neg_start && pos_start == neg_end && pos_end > neg_end) {
        fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
        ovp = false;
        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = get_cov_sign(0, pos_start-refstart, pos_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, pos_start[%d]-refstart, pos_end[%d]-refstart, bpcov) (%f) = %f\n",pos_start, pos_end, get_cov_sign(0, pos_start-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = 0;   
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);         
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            if (chaining_hold_strand == '-' && last_ovp_end > pos_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, neg_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov (%f) + get_cov_sign(0, last_ovp_end[%d]-refstart, neg_end[%d]-refstart, bpcov) (%f) = %f\n", chaining_hold_cov, last_ovp_end, neg_end, get_cov_sign(0, last_ovp_end-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, neg_end[%d]-refstart, bpcov) (%f) = %f\n", neg_start, neg_end, get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            }
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            last_ovp_end = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }
    } else if (pos_start > neg_end) {
        fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
        ovp = false;
        if (target_node_strand == '+') {
            //Positive gonna be pushed to next node
            end_chaining_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
            fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, pos_start[%d]-refstart, pos_end[%d]-refstart, bpcov) (%f) = %f\n",pos_start, pos_end, get_cov_sign(0, pos_start-refstart, pos_end-refstart, bpcov), end_chaining_cov);
            chaining_hold_strand = '-';
            chaining_hold_cov = 0;
            last_ovp_end = 0;       
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);     
        } else if (target_node_strand == '-') {
            //Negative gonna be pushed to next node
            if (chaining_hold_strand == '-' && last_ovp_end > neg_start) {
                end_chaining_cov = chaining_hold_cov + get_cov_sign(0, last_ovp_end-refstart, neg_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: chaining_hold_cov (%f) + get_cov_sign(0, last_ovp_end[%d]-refstart, neg_end[%d]-refstart, bpcov) (%f) = %f\n", chaining_hold_cov, last_ovp_end, neg_end, get_cov_sign(0, last_ovp_end-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            } else {
                end_chaining_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
                fprintf(stderr, "\t\t** end_chaining_cov: get_cov_sign(0, neg_start[%d]-refstart, neg_end[%d]-refstart, bpcov) (%f) = %f\n", neg_start, neg_end, get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov), end_chaining_cov);
            }
            chaining_hold_strand = '+';
            chaining_hold_cov = 0;
            last_ovp_end = 0;
            fprintf(stderr, "\t\t** chaining_hold_cov: %f \n", chaining_hold_cov);
        }
    }

    fprintf(stderr, ">>> ovp : %f \n", ovp);
    fprintf(stderr, ">>> Final coverage (%c): %f \n", target_node_strand, end_chaining_cov);
    fprintf(stderr, ">>> chaining_hold_cov (%c): %f, last_ovp_end: %d \n", chaining_hold_strand, chaining_hold_cov, last_ovp_end);

    if (end_chaining_cov < 0.5) end_chaining_cov=0.0;
    if (chaining_hold_cov < 0.5) chaining_hold_cov = 0.0;

    return ovp;
}




void ovp_coverage_push_node(int& g_idx, int& g_num, int& n_idx, int& n_num, bool& reach_end) {
    if (g_idx < g_num-1) {
        if (n_idx < n_num-1) {
            n_idx++;
        } else if (n_idx == n_num-1) {
            g_idx += 1;
            n_idx = 1;
        }
    } else if (g_idx == g_num-1) {
        if (n_idx < n_num-1) {
            n_idx++;
        } else if (n_idx == n_num-1) {
            reach_end = true;
        }
    }
}


void calculate_ovp_coverage(int pos_start, int pos_end, int neg_start, int neg_end, int refstart, int refend, float bundle_coverage_ratio_pos, float bundle_coverage_ratio_neg, float& pos_cov, float& neg_cov, GVec<float>* bpcov) {
    // fprintf(stderr, "Inside `segs_overlap`\n");
    fprintf(stderr, "Process read (pre: %d - %d ;  now: %d - %d)!!!\n", pos_start, pos_end, neg_start, neg_end);
    // fprintf(stderr, "\t>>>>> brec: %d - %d\n", pos_start, pos_end);

    // node_coverage_neg = get_cov(0, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
    // node_coverage_uns = get_cov(1, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
    // node_coverage_pos = get_cov(2, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
    float pos_cov_ovp = 0.0;
    float uns_cov_ovp = 0.0;
    float neg_cov_ovp = 0.0;
    if (pos_end <  neg_start) {
        // ----------   |(s).................(e)|
        // The read is outside the current bundle => skipped!
        fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
        pos_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
    } else if (pos_start < neg_start && pos_end == neg_start) {
        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
        fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
        pos_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
    } else if (pos_start < neg_start && pos_end > neg_start && pos_end < neg_end) {
        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
        fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");

        pos_cov_ovp = get_cov(2, neg_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, pos_end-refstart, bpcov);

        pos_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov);
        neg_cov = get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
    } else if (pos_start < neg_start && pos_end > neg_start && pos_end == neg_end) {
        // -----|(s)---------(e)|
        fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
        pos_cov_ovp = get_cov(2, neg_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, pos_end-refstart, bpcov);

        pos_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov);
        neg_cov = 0;
        // get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
    } else if (pos_start < neg_start && pos_end > neg_start && pos_end > neg_end) {
        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
        fprintf(stderr, "\t** Bundle: -----|(s)----------(e)|------\n");

        pos_cov_ovp = get_cov(2, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);

        pos_cov = get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov);
        pos_cov += get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov);
        neg_cov = 0;
        // get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
    } else if (pos_start == neg_start && pos_start < neg_end && pos_end < neg_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
        pos_cov_ovp = get_cov(2, neg_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, pos_end-refstart, bpcov);

        pos_cov = 0;
        // get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov);
        neg_cov = get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);

    } else if (pos_start == neg_start && pos_start < neg_end && pos_end == neg_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
        pos_cov_ovp = get_cov(2, neg_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, pos_end-refstart, bpcov);

        pos_cov = 0;
        // get_cov_sign(2, pos_start-refstart, neg_start-refstart, bpcov);
        neg_cov = 0;
        // get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);

    } else if (pos_start == neg_start && pos_start < neg_end && pos_end > neg_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
        pos_cov_ovp = get_cov(2, neg_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, neg_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);

        pos_cov = get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov);
        neg_cov = 0;
        // get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov);

    } else if (pos_start > neg_start && pos_start < neg_end && pos_end > neg_start && pos_end < neg_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");

        fprintf(stderr, "\t** neg_start: %d;  pos_start: %d;  pos_end: %d;  neg_end: %d;\n", neg_start, pos_start, pos_end, neg_end);


        pos_cov_ovp = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);

        pos_cov = 0;
        // get_cov_sign(2, pos_start-refstart, pos_start-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov);

        fprintf(stderr, "\t** Just checking@!!! neg_cov: %f\n", neg_cov);
        neg_cov += get_cov_sign(0, pos_end-refstart, neg_end-refstart, bpcov);
        fprintf(stderr, "\t** Just checking@!!! neg_cov: %f\n", neg_cov);



        float test_neg_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
        float test_neg_cov_only = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);

        fprintf(stderr, "\t** Whole region!!! neg_cov: %f;  test_neg_cov_only: %f\n", test_neg_cov, test_neg_cov_only);



    } else if (pos_start > neg_start && pos_start < neg_end && pos_end == neg_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");

        fprintf(stderr, "\t** neg_start: %d;  pos_start: %d;  pos_end: %d;  neg_end: %d;\n", neg_start, pos_start, pos_end, neg_end);

        pos_cov_ovp = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);

        pos_cov = 0;
        // get_cov_sign(2, pos_start-refstart, pos_start-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov);
        fprintf(stderr, "\t** Just checking@!!! neg_cov: %f\n", neg_cov);

        float test_neg_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
        float test_neg_cov_only = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);

        fprintf(stderr, "\t** Whole region!!! neg_cov: %f;  test_neg_cov_only: %f\n", test_neg_cov, test_neg_cov_only);


    } else if (pos_start > neg_start && pos_start < neg_end && pos_end > neg_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
        pos_cov_ovp = get_cov(2, pos_start-refstart, neg_end-refstart, bpcov);
        uns_cov_ovp = get_cov(1, pos_start-refstart, neg_end-refstart, bpcov);
        neg_cov_ovp = get_cov(0, pos_start-refstart, neg_end-refstart, bpcov);

        pos_cov = get_cov_sign(2, neg_end-refstart, pos_end-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, pos_start-refstart, bpcov);
    } else if (pos_start > neg_start && pos_start == neg_end && pos_end > neg_end) {
        // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
        // The overlapping with the current processing bundle.
        fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
        pos_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);

    } else if (pos_start > neg_end) {
        fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
        pos_cov = get_cov_sign(2, pos_start-refstart, pos_end-refstart, bpcov);
        neg_cov = get_cov_sign(0, neg_start-refstart, neg_end-refstart, bpcov);
    }

    uns_cov_ovp = uns_cov_ovp - pos_cov_ovp - neg_cov_ovp;
    
    fprintf(stderr, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
    fprintf(stderr, "&&&& Starting cov redistribution &&&&&&&&&&&&&&\n");
    fprintf(stderr, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");

    // float pos_cov_pos = get_cov(2, pos_start-refstart, pos_end-refstart, bpcov);
    // float pos_cov_uns = get_cov(1, pos_start-refstart, pos_end-refstart, bpcov);
    // float pos_cov_neg = get_cov(0, pos_start-refstart, pos_end-refstart, bpcov);
    // pos_cov_uns = pos_cov_uns - pos_cov_pos - pos_cov_neg;

    // float neg_cov_pos = get_cov(2, neg_start-refstart, neg_end-refstart, bpcov);
    // float neg_cov_uns = get_cov(1, neg_start-refstart, neg_end-refstart, bpcov);
    // float neg_cov_neg = get_cov(0, neg_start-refstart, neg_end-refstart, bpcov);
    // neg_cov_uns = neg_cov_uns - neg_cov_pos - neg_cov_neg;

    // fprintf(stderr, "pos_cov_pos: %f, pos_cov_uns: %f, pos_cov_neg: %f\n", pos_cov_pos, pos_cov_uns, pos_cov_neg);
    // fprintf(stderr, "neg_cov_pos: %f, neg_cov_uns: %f, neg_cov_neg: %f\n", neg_cov_pos, neg_cov_uns, neg_cov_neg);

    fprintf(stderr, "pos_cov_ovp: %f, neg_cov_ovp: %f, uns_cov_ovp: %f\n", pos_cov_ovp, neg_cov_ovp, uns_cov_ovp);

    fprintf(stderr, "++++++ positive decomposition: pos_cov: %f, pos_cov_ovp: %f, uns_cov_ovp*bundle_coverage_ratio_pos: %f\n", pos_cov, pos_cov_ovp, uns_cov_ovp*bundle_coverage_ratio_pos);
    fprintf(stderr, "------ negative decomposition: neg_cov: %f, neg_cov_ovp: %f, uns_cov_ovp*bundle_coverage_ratio_neg: %f\n", neg_cov, neg_cov_ovp, uns_cov_ovp*bundle_coverage_ratio_neg);


    pos_cov = pos_cov + pos_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_pos;
    neg_cov = neg_cov + neg_cov_ovp + uns_cov_ovp*bundle_coverage_ratio_neg;

    // fprintf(stderr, "pos_ratio: %f, neg_ratio: %f, pos_dis_cov: %f, neg_dis_cov: %f \n", 0.5, 0.5, uns_cov_ovp/2, uns_cov_ovp/2);



    // if (pos_cov_ovp==0 && neg_cov_ovp==0) {
    //     pos_cov = pos_cov + pos_cov_ovp + uns_cov_ovp/2;
    //     neg_cov = neg_cov + neg_cov_ovp + uns_cov_ovp/2;
    //     fprintf(stderr, "pos_ratio: %f, neg_ratio: %f, pos_dis_cov: %f, neg_dis_cov: %f \n", 0.5, 0.5, uns_cov_ovp/2, uns_cov_ovp/2);
    // } else {
    //     pos_cov = pos_cov + pos_cov_ovp + uns_cov_ovp*(pos_cov_ovp/(pos_cov_ovp+neg_cov_ovp));
    //     neg_cov = neg_cov + neg_cov_ovp + uns_cov_ovp*(neg_cov_ovp/(pos_cov_ovp+neg_cov_ovp));
    //     fprintf(stderr, "pos_ratio: %f, neg_ratio: %f, pos_dis_cov: %f, neg_dis_cov: %f \n", pos_cov_ovp/(pos_cov_ovp+neg_cov_ovp), neg_cov_ovp/(pos_cov_ovp+neg_cov_ovp), uns_cov_ovp*(pos_cov_ovp/(pos_cov_ovp+neg_cov_ovp)), uns_cov_ovp*(neg_cov_ovp/(pos_cov_ovp+neg_cov_ovp)));

    // }
    fprintf(stderr, "\n");
}

void redistribute_unstranded_rcov(float* rprop, GVec<float>* bpcov, int refstart, int refend, int rstart, int rend) {

    float read_coverage_neg = get_cov(0, rstart-refstart+1, rend-refstart+1, bpcov);
	float read_coverage_uns = get_cov(1, rstart-refstart+1, rend-refstart+1, bpcov);
	float read_coverage_pos = get_cov(2, rstart-refstart+1, rend-refstart+1, bpcov);
	read_coverage_uns = read_coverage_uns - read_coverage_neg - read_coverage_pos;

    if (read_coverage_neg == 0.0 && read_coverage_pos == 0.0) {
        rprop[0] = 0.5;
        rprop[1] = 0.5;
    } else if (read_coverage_neg == 0.0 && read_coverage_pos != 0.0) {
        rprop[0] = 0.0;
        rprop[1] = 1.0;
    } else if (read_coverage_neg != 0.0 && read_coverage_pos == 0.0) {
        rprop[0] = 1.0;
        rprop[1] = 0.0;
    } else {
        rprop[0] = read_coverage_neg/(read_coverage_neg+read_coverage_pos);
        rprop[1] = 1.0 - rprop[0];
    }

	fprintf(stderr, ">> Ref  (%d - %d)\n", refstart, refend);
	fprintf(stderr, ">> read (%u - %u)\n", rstart, rend);
    if (read_coverage_neg != 0.0 && read_coverage_pos != 0.0) {
        fprintf(stderr, ">> read_coverage_neg: %f \n", read_coverage_neg);

        fprintf(stderr, ">> read_coverage_uns: %f \n", read_coverage_uns);

        fprintf(stderr, ">> read_coverage_pos: %f \n", read_coverage_pos);         
    }


	// if(readlist[n]->strand != 1 && readlist[n]->strand != -1) {
	// 	fprintf(stderr, "UNSTRANDED!!!\n");
	// }
}

int calOverlapLen(uint rstart, uint rend, uint start, uint end) {
    fprintf(stderr, ">> Inside calOverlapLen:!!!\n");
    fprintf(stderr, ">> (%u - %u); (%u - %u)!!!\n", rstart, rend, start, end);

    if (rstart>rend) { Gswap(rstart,rend); }
    if (start<rstart) {
        if (rstart>end) return 0;
        return (rend>end) ? end-rstart+1 : rend-rstart+1;
    }
    else { //rstart<=start
        if (start>rend) return 0;
        return (rend<end)? rend-start+1 : end-start+1;
    }
}

char* sprintTime() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}

int loadPtFeatures(FILE* f, GArray<GRefPtData>& refpts) {
  //expected format:
  //<chromosome> <coordinate> <strand> <feature_type>
  int num=0;
  GLineReader lr(f);
  char* line=NULL;
  GDynArray<char*> tokens;
  while ((line=lr.nextLine())!=NULL) {
    strsplit(line, tokens);
    if (tokens.Count()<4)
    	GError("Error parsing point-feature line (not enough columns):\n%s\n",line);
    int start;
    if (!strToInt(tokens[1], start))
    	GError("Error parsing point-feature line (invalid coordinate):\n%s\n",line);
    int8_t strand=-2;
    if (strlen(tokens[2])==1) {
    	if (tokens[2][0]=='+')
    		strand=1;
    	else if (tokens[2][0]=='-')
    		strand=-1;
    	else if (tokens[2][0]=='.')
    		strand=0;
    }
    if (strand==-2)
    	GError("Error parsing point-feature line (invalid strand):\n%s\n",line);
    GPFType pftype=GPFT_NONE;
    if (strcmp(tokens[3], "TSS")==0)
			pftype=GPFT_TSS;
	else if (strcmp(tokens[3], "CPAS")==0)
			pftype=GPFT_CPAS;
    if (pftype==0)
    	GError("Error parsing point-feature line (unrecognized type):\n%s\n",line);
    GPtFeature* ptf=new GPtFeature(pftype, -1, strand, start);
    addPtFeature(tokens[0], ptf, refpts);
    num++;
  } //while line
  return num;
}


void addPtFeature(const char* refname, GPtFeature* pf, GArray<GRefPtData>& refpts) {
  //expects gseqNames to be set to GffObj::names and initialized/populated already!
  //MUST be called AFTER the guides file has been loaded (if given)
  int gseq_id=gseqNames->gseqs.addName(refname);
  pf->ref_id=gseq_id;
  int ridx=-1;
  GRefPtData* rpd=NULL;
  GRefPtData rd(gseq_id);
  if (refpts.Count()==0) {
	  ridx=refpts.Add(rd);
  } else {
	  ridx=refpts.IndexOf(rd);
	  if (ridx<0) {
		  ridx=refpts.Add(rd);
	  }
  }
  if (ridx<0) GError("Error adding GRefPtData entry (bug!)\n");
  rpd = & refpts.Get(ridx);
  rpd->add(pf);
}
