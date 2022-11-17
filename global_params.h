#ifndef __GLOBAL_PARAMS_H__
#define __GLOBAL_PARAMS_H__

#pragma once

#include "definitions.h"
// #include "dot_record.h"
#include "mode.hpp"
// #include "rlink.h"
// #include "t_record.h"
// #include "tablemaker.h"
// #include "unispg.h"

#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "GThreads.h"
#include "GHashMap.hh"
#include "GFaSeqGet.h"
#include "GFaSeqGet.h"
#include "gff.h"
#include "GSam.h"

#include <fstream>
#include <iostream>
using namespace std;
/*******************************************
 ** Argument parsing parameters.
 *******************************************/

// class GFastaDb;
// class GStr;
// template <typename K> struct GHashKey_wyHash;
// template <typename K> struct GHashKey_Eq;
// template <class Hash=GHashKey_wyHash<const char*>, class Eq=GHashKey_Eq<const char*>, typename khInt_t=uint64_t> class GStrSet;
// struct FILE;
struct TInputFiles;
struct DOTInputFile;
// class GffNames;

extern bool debugMode; // "debug" or "D" tag.
extern bool verbose; // "verbose" / "v" tag.
extern bool multiMode; // "multi" tag.
extern bool graph_bed;
extern bool keepTempFiles; // "keeptmp" tag.
extern GFastaDb* gfasta; // "rseq" or "S" tag.
extern GStr ptff; // "ptf" tag. (point features)
extern bool fr_strand; // "fr" tag.
extern bool rf_strand; // "rf" tag.
extern bool includesource; // "z" tag.
extern bool retained_intron; // "i" tag. set by parameter -i for merge option
extern bool trim; // "t" tag. 
extern bool eonly; // "e" tag. for mergeMode includes estimated coverage sum in the merged transcripts
extern bool nomulti; // "u" tag.
extern GStrSet<> excludeGseqs; // "x" tag. hash of chromosomes/contigs to exclude (e.g. chrM)
extern bool guided; // "G" tag.
extern GStr guidegff; // "G" tag
extern GStr label; // "l" tag.
extern int mintranscriptlen; // "m" tag. minimum length for a transcript to be printed
extern uint junctionsupport; // "a" tag. anchor length for junction to be considered well supported <- consider shorter??
extern int junctionthr; // "j" tag. number of reads needed to support a particular junction
extern float singlethr; // "s" tag. coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
extern float readthr; // "c" tag. read coverage per bundle bp to accept it; // paper uses 3
extern float isofrac; // "f" tag. 
extern int num_cpus; // "p" tag.
extern uint bundledist;  // "g" tag. reads at what distance should be considered part of separate bundles
extern uint runoffdist; // threshold for 'bundledist'
extern float mcov; // "M" tag. fraction of bundle allowed to be covered by multi-hit reads paper uses 1
extern bool geneabundance; // "A" tag.
extern uint sserror; // "E" tag. window arround splice sites that we use to generate consensus in case of long read data
extern float fpkm_thr; // "F" tag.
extern float tpm_thr; // "T" tag.
extern bool isunitig; // "U" tag.
extern bool enableNames; // "E" tag (removed)
extern FILE* c_out; // "C" tag.
extern FILE* f_out; // get from "outfname" param.


/*******************************************
 ** File-related parameters.
 *******************************************/
// Dot file => for reading in DOT format.
extern GStr unispgdotfname_root; 
extern GStr unispgdotfname;
extern GStr unispgdotfname_pos; 
extern GStr unispgdotfname_neg; 
// Annotation file => outputing in GFF format.
extern GStr outfname;
extern GStr f_basename;
extern GStr outfname_prefix;
extern GStr out_dir;
extern GStr tmp_path;
extern GStr cram_ref; //"ref" / "cram-ref" tag. Reference genome FASTA for CRAM input
extern GStr tmpfname; // "o" tag.
extern GStr genefname;
// Ratio file => for coverage comparison & visualization.
extern ofstream cov_file_pos;
extern ofstream cov_file_neg;
extern ofstream cov_file_pos_norm;
extern ofstream cov_file_neg_norm;

/*******************************************
 ** Program-defined global parameter.
 *******************************************/
extern multiStringTieMode mode;
extern bool NoMoreBundles;
extern GffNames* gseqNames; //used as a dictionary for reference sequence names and ids
extern int refseqCount; // number of reference sequences found in the guides file
extern bool includecov;
// Global counter variable.
extern int sample_num;
extern double Num_Fragments; //global fragment counter (aligned pairs)
extern double Frag_Len;
extern double Cov_Sum;
extern int GeneNo; //-- global "gene" counter
// For multistringtie applyUNISPG only. Counting how many reads are skipped (not processed).
extern int skip_counter;
extern int boundary_counter;
extern int skip_counter_nh;
extern int allowed_nodes;

/*******************************************
 ** Reader parameters.
 *******************************************/
extern TInputFiles bamreader;
extern DOTInputFile dotreader;
extern DOTInputFile dotreader_pos;
extern DOTInputFile dotreader_neg;



#ifndef NOTHREADS
//single producer, multiple consumers
//main thread/program is always loading the producer
extern GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
extern GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
extern GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
extern int bundleWork; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue

//GFastMutex waitMutex;
extern GMutex waitMutex; // controls threadsWaiting (idle threads counter)

extern int threadsWaiting; // idle worker threads
extern GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle

extern GConditionVar haveClear; //will notify when bundle buf space available

extern GMutex queueMutex; //controls bundleQueue and bundles access

extern GFastMutex printMutex; //for writing the output to file

extern GFastMutex logMutex; //only when verbose - to avoid mangling the log output

extern GFastMutex bamReadingMutex;

extern GFastMutex countMutex;
#endif

#endif