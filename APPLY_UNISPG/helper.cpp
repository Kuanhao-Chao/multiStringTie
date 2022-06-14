#include "helper.h"

bool segs_overlap(int s1_start, int s1_end, int s2_start, int s2_end) {
    bool ovp =true;
    // fprintf(stderr, "Process read (pre: %d - %d ;  now: %d - %d)!!!\n", pre_refstart, pre_refend, s1_start, s1_end);
    // fprintf(stderr, "\t>>>>> brec: %d - %d\n", s1_start, s1_end);
    if (s1_end <  s2_start) {
        // ----------   |(s).................(e)|
        // The read is outside the current bundle => skipped!
        fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
        ovp = false;
    } else if (s1_start < s2_start && s1_end == s2_start) {
        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
        fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
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
    } else if (s1_start > s2_start && s1_start < s2_end && s1_end > s2_end) {
        // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
        // The overlapping with the current processing bundle.
        fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
    } else if (s1_start > s2_start && s1_start == s2_end && s1_end > s2_end) {
        // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
        // The overlapping with the current processing bundle.
        fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
    } else if (s1_start > s2_end) {
        fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
        ovp = false;
    }

    return ovp;
}