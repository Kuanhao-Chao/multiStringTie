#pragma once
#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "rlink.h"
#include "t_record.h"
#include "dot_record.h"
#include "unispg.h"
#include "helper.h"
#include "definitions.h"

#include "multistringtie_APPLY.cpp"
#include "multistringtie_CREATE.cpp"
// Include fstream to quickly write out the ratio!!!
#include <iostream>
#include <fstream>
#include <vector> 
#include "time.h"

/*******************************************
 ** Argument parsing parameters.
 *******************************************/
bool debugMode=false; // "debug" or "D" tag.
bool verbose=false; // "verbose" / "v" tag.
bool ballgown=false; // "B" tag.
GStr ballgown_dir; // "b" tag.
bool viral=false; // "viral" tag.
bool mixedMode=false; // "mix" tag. both short and long read data alignments are provided
bool mergeMode = false; // "merge" tag. For running StringTie Merge.
bool multiMode=false; // "multi" tag.
bool graph_bed=false; // "graph_bed" tag.
bool keepTempFiles; // "keeptmp" tag.
GFastaDb* gfasta=NULL; // "rseq" or "S" tag.
GStr ptff; // "ptf" tag. (point features)
bool fr_strand=false; // "fr" tag.
bool rf_strand=false; // "rf" tag.
bool includesource=true; // "z" tag.
bool retained_intron=false; // "i" tag. set by parameter -i for merge option
bool trim=true; // "t" tag. 
bool eonly=false; // "e" tag. for mergeMode includes estimated coverage sum in the merged transcripts
bool nomulti=false; // "u" tag.
bool longreads=false; // "L" tag.
bool rawreads=false; // "R" tag.
GStrSet<> excludeGseqs; // "x" tag. hash of chromosomes/contigs to exclude (e.g. chrM)
bool guided=false; // "G" tag.
GStr guidegff; // "G" tag
GStr label("STRG"); // "l" tag.
int mintranscriptlen=200; // "m" tag. minimum length for a transcript to be printed
uint junctionsupport=10; // "a" tag. anchor length for junction to be considered well supported <- consider shorter??
int junctionthr=1; // "j" tag. number of reads needed to support a particular junction
float singlethr=4.75; // "s" tag. coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
float readthr=1; // "c" tag. read coverage per bundle bp to accept it; // paper uses 3
float isofrac=0.01; // "f" tag. 
int num_cpus=1; // "p" tag.
uint bundledist=50;  // "g" tag. reads at what distance should be considered part of separate bundles
uint runoffdist=200; // threshold for 'bundledist'
float mcov=1; // "M" tag. fraction of bundle allowed to be covered by multi-hit reads paper uses 1
bool geneabundance=false; // "A" tag.
uint sserror=25; // "E" tag. window arround splice sites that we use to generate consensus in case of long read data
float fpkm_thr=1; // "F" tag.
float tpm_thr=1; // "T" tag.
bool isunitig=true; // "U" tag.
bool enableNames=false; // "E" tag (removed)
FILE* c_out=NULL; // "C" tag.
FILE* f_out=NULL; // get from "outfname" param.


/*******************************************
 ** File-related parameters.
 *******************************************/
// Dot file => for reading in DOT format.
GStr unispgdotfname_root; 
GStr unispgdotfname_pos; 
GStr unispgdotfname_neg; 
// Annotation file => outputing in GFF format.
GStr outfname;
GStr outfname_prefix;
GStr out_dir;
GStr f_basename;
GStr tmp_path;
GStr cram_ref; //"ref" / "cram-ref" tag. Reference genome FASTA for CRAM input
GStr tmpfname; // "o" tag.
GStr genefname;
GStr traindir; // "cds" tag. training directory for CDS option (removed)
// Ratio file => for coverage comparison & visualization.
ofstream cov_file_pos;
ofstream cov_file_neg;
ofstream cov_file_pos_norm;
ofstream cov_file_neg_norm;

/*******************************************
 ** Program-defined global parameter.
 *******************************************/
multiStringTieMode mode;
bool NoMoreBundles=false;
GffNames* gseqNames=NULL; //used as a dictionary for reference sequence names and ids
int refseqCount=0; // number of reference sequences found in the guides file
bool includecov=false;
// Global counter variable.
int sample_num = 1;
double Num_Fragments=0; //global fragment counter (aligned pairs)
double Frag_Len=0;
double Cov_Sum=0;
int GeneNo=0; //-- global "gene" counter
// For multistringtie applyUNISPG only. Counting how many reads are skipped (not processed).
int skip_counter = 0;
int boundary_counter = 0;
int skip_counter_nh = 0;
int allowed_nodes=1000;

/*******************************************
 ** Reader parameters.
 *******************************************/
TInputFiles bamreader;
DOTInputFile dotreader;
DOTInputFile dotreader_pos;
DOTInputFile dotreader_neg;

/*******************************************
 ** Multithread parameters.
 *******************************************/
#ifndef NOTHREADS
//single producer, multiple consumers
//main thread/program is always loading the producer
GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
int bundleWork=1; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue
//GFastMutex waitMutex;
GMutex waitMutex; // controls threadsWaiting (idle threads counter)
int threadsWaiting; // idle worker threads
GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle
GConditionVar haveClear; //will notify when bundle buf space available
GMutex queueMutex; //controls bundleQueue and bundles access
GFastMutex printMutex; //for writing the output to file
GFastMutex logMutex; //only when verbose - to avoid mangling the log output
GFastMutex bamReadingMutex;
GFastMutex countMutex;
#endif

int main(int argc, char*argv[]) {
	std::cout << "This is the multiStringTie program.\n" << std::endl;

	fprintf(stderr, "%d  %s  %s  %s  %s\n", argc, argv[0], argv[1], argv[2], argv[3]);
	if (strcmp(argv[1], "CREATE_UNISPG") || strcmp(argv[1], "APPLY_UNISPG")) {
		// fprintf(stderr, "This is mode: %s  %d\n", argv[1], strcmp(argv[1], "APPLY_UNISPG"));
		if (strcmp(argv[1], "CREATE_UNISPG") == 0) {
			mode = CREATE_UNISPG;
		}
		if (strcmp(argv[1], "APPLY_UNISPG") == 0) {
			mode = APPLY_UNISPG;
		}
		for (int i=1; i<argc-1; i++) {
			argv[i] = argv[i+1];
		}
		argv[argc-1] = NULL;
		argc = argc-1;
		fprintf(stderr, "%d  %s  %s  %s  %s\n", argc, argv[0], argv[1], argv[2], argv[3]);
	}
   
	if (mode == CREATE_UNISPG) {
		multistringtie_CREATE(argc, argv);
	} else if (mode == APPLY_UNISPG) {
		multistringtie_APPLY(argc, argv);
	}
    return 0;
}