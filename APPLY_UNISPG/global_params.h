#include "mode.hpp"
#include "GHashMap.hh"
#include "GFaSeqGet.h"
// #include "tmerge.h"

/*******************************************
 ** Argument parsing parameters.
 *******************************************/
extern bool debugMode; // "debug" or "D" tag.
extern bool verbose; // "v" tag.
extern bool ballgown; // "b" tag.
extern bool keepTempFiles; // "keeptmp" tag.
extern bool fr_strand; // "fr" tag.
extern bool rf_strand; // "rf" tag.
extern bool guided; // "G" tag.
extern bool viral; // "viral" tag.
extern bool eonly; // "e" tag. for mergeMode includes estimated coverage sum in the merged transcripts
extern bool longreads; // "L" tag.
extern GStrSet<> excludeGseqs; // "x" tag. hash of chromosomes/contigs to exclude (e.g. chrM)
extern uint bundledist;  // "g" tag. reads at what distance should be considered part of separate bundles
extern uint runoffdist; // threshold for 'bundledist'
extern float readthr; // "c" tag. read coverage per bundle bp to accept it; // paper uses 3
extern float tpm_thr; // "T" tag.
extern float fpkm_thr; // "F" tag.
extern bool isunitig; // "U" tag.
extern uint junctionsupport; // "a" tag. anchor length for junction to be considered well supported <- consider shorter??
extern uint sserror; // "E" tag. window arround splice sites that we use to generate consensus in case of long read data
extern int junctionthr; // "j" tag. number of reads needed to support a particular junction
extern float mcov; // "M" tag. fraction of bundle allowed to be covered by multi-hit reads paper uses 1
extern int mintranscriptlen; // "m" tag. minimum length for a transcript to be printed
extern float singlethr; // "s" tag. coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
extern bool mixedMode; // "mix" tag. both short and long read data alignments are provided
extern float isofrac; // "f" tag. 
extern bool trim; // "t" tag. 
extern bool includesource; // "z" tag.
extern bool nomulti; // "u" tag.
extern bool multiMode; // "multi" tag.
extern bool mergeMode; // "merge" tag. 
extern GFastaDb* gfasta; // "rseq" or "S" tag.
extern int num_cpus; // "p" tag.
extern bool rawreads; // "R" tag.
extern GStr label; // "l" tag.
extern bool retained_intron; // "i" tag. set by parameter -i for merge option
extern bool geneabundance; // "A" tag.
extern GStr guidegff; // "G" tag
extern GStr ptff; // "ptf" tag. (point features)
extern FILE* f_out; // get from "outfname" param.
extern FILE* c_out; // "C" tag.

/*******************************************
 ** File-related parameters.
 *******************************************/
// Dot file => for reading in DOT format.
extern GStr unispgdotfname_root; 
extern GStr unispgdotfname_pos; 
extern GStr unispgdotfname_neg; 
// Annotation file => outputing in GFF format.
extern GStr outfname;
extern GStr out_dir;
extern GStr tmp_path;
extern GStr cram_ref; //reference genome FASTA for CRAM input
extern GStr tmpfname;
extern GStr genefname;
extern GStr traindir; // training directory for CDS option
// Ratio file => for coverage comparison & visualization.
// extern ofstream ratio_file_pos;
// extern ofstream ratio_file_neg;

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

/*******************************************
 ** Reader parameters.
 *******************************************/
// extern TInputFiles bamreader;
// extern DOTInputFile dotreader;
// extern DOTInputFile dotreader_pos;
// extern DOTInputFile dotreader_neg;