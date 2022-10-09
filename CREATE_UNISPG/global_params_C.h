#ifndef __GLOBAL_PARAMS_C_H__
#define __GLOBAL_PARAMS_C_H__

#pragma once

// #include "unispg_C.h"
#include "../global_params.h"
#include "../rlink.h"
// struct GRefData;
// struct GRefPtData;
// struct RC_TData;
// struct RC_Feature;
// class UnispgGp_CREATE;

/*******************************************
 ** CREATE_UNISPG specific parameters.
 *******************************************/
// extern UnispgGp_CREATE* unispg_gp;
extern bool universal_splice_graph;
extern FILE* uinigraph_out;
extern GStr unigraphfname; 
extern GStr plot_dir;

/*****************************
 * Declaring DOT file.
 *****************************/
extern GVec<FILE*> pos_dot_vec;
extern GVec<GStr> pos_dotfname_vec; 
extern GVec<FILE*> neg_dot_vec;
extern GVec<GStr> neg_dotfname_vec; 
extern FILE* dot;
extern GStr dotfname; 
extern GVec<FILE*>* dot_vec[2];
extern GVec<GStr>* dotfname_vec[2]; 

/*****************************
 * Declaring BED file.
 *****************************/
extern GVec<FILE*>* node_lclg_bed_vec[2];
extern GVec<GStr>* nodelclgfname_vec[2]; 
extern GVec<FILE*>* edge_lclg_bed_vec[2];
extern GVec<GStr>* edgelclgfname_vec[2]; 

extern GVec<FILE*>* node_novp_bed_vec[2];
extern GVec<GStr>* nodenovpfname_vec[2]; 
extern GVec<FILE*>* edge_novp_bed_vec[2];
extern GVec<GStr>* edgenovpfname_vec[2]; 

extern GVec<FILE*>* node_unispg_bed_vec[2];
extern GVec<GStr>* nodeunispgfname_vec[2]; 
extern GVec<FILE*>* edge_unispg_bed_vec[2];
extern GVec<GStr>* edgeunispgfname_vec[2]; 

extern FILE* node_unispg_unstrand_bed; 
extern GStr nodeunispgfname_unstrand; 
// These are for temporary file holding
extern FILE* node_cov_bed;
extern GStr nodecovfname; 
extern FILE* edge_cov_bed;
extern GStr edgecovfname; 

/*****************************
 * Declaring reference related data structure
 *****************************/
extern GVec<GRefData> refguides; // plain vector with transcripts for each chromosome
extern GArray<GRefPtData> refpts; // sorted,unique array of refseq point-features data

/*****************************
 * Declaring Ballgown related data structure
 *   table indexes for Ballgown Raw Counts data (-B/-b option)
 *****************************/
extern GPVec<RC_TData> guides_RC_tdatarlink_C; //raw count data or other info for all guide transcripts
extern GPVec<RC_Feature> guides_RC_exonsrlink_C; //raw count data for all guide exons
extern GPVec<RC_Feature> guides_RC_intronsrlink_C;//raw count data for all guide introns

#endif