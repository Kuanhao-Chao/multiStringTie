#include "helper.h"

bool segs_overlap(int s1_start, int s1_end, int s2_start, int s2_end) {
    fprintf(stderr, "Inside `segs_overlap`\n");
    bool ovp =true;
    fprintf(stderr, "Process two segments (pre: %d - %d ;  now: %d - %d)!!!\n", s1_start, s1_end, s2_start, s2_end);
    // fprintf(stderr, "\t>>>>> brec: %d - %d\n", s1_start, s1_end);
    if (s1_end <  s2_start) {
        // ----------   |(s).................(e)|
        // The read is outside the current bundle => skipped!
        fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
        ovp = false;
    } else if (s1_start < s2_start && s1_end == s2_start) {
        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
        fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
        // ovp = false;
    } else if (s1_start < s2_start && s1_end > s2_start) {
        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
        fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
    } else if (s1_start < s2_start && s1_end == s2_start) {
        // -----|(s)---------(e)|
        fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
    } else if (s1_start == s2_start && s1_start < s2_end && s1_end < s2_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
    } else if (s1_start == s2_start && s1_start < s2_end && s1_end == s2_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
    } else if (s1_start == s2_start && s1_start < s2_end && s1_end > s2_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end < s2_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end == s2_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end > s2_end) {
        // |(s)----------.................(e)|   or   |(s)....----------........(e)|
        fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
    } else if (s1_start > s2_start && s1_start == s2_end && s1_end > s2_end) {
        // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
        // The overlapping with the current processing bundle.
        fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
        // ovp = false;
    } else if (s1_start > s2_end) {
        fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
        ovp = false;
    }

    return ovp;
}


void calculate_ovp_coverage(int pos_start, int pos_end, int neg_start, int neg_end, int refstart, int refend, float bundle_coverage_ratio_pos, float bundle_coverage_ratio_neg, float& pos_cov, float& neg_cov, GVec<float>* bpcov) {
    fprintf(stderr, "Inside `segs_overlap`\n");
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
    } else if (pos_start < neg_start && pos_end > neg_start && pos_end == neg_start) {
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