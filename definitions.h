#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__

// Let's test single thread first!! 
#define NOTHREADS

//#define GMEMTRACE 1 
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define VERSION "0.1.0"

#define USAGE "multiStringTie v" VERSION " usage:\n\n\
multiStringTie [CREATE_UNISPG / APPLY_UNISPG] <in.DOT ..> <in.bam ..> \
 [--mix] [--conservative] [--rf] [--fr]\n\
Assemble RNA-Seq alignments into potential transcripts with universal splice graph.\n\
Options:\n\
* Transcript CREATE_UNISPG mode: \n\
   multistringtie CREATE_UNISPG --multi [Options] { bam_list | sample1.bam ...}\n\
    --version : print just the version at stdout and exit\n\
    -o output path/file name for the assembled transcripts GTF (default: stdout)\n\n\
* APPLY_UNISPG mode: \n\
   multistringtie APPLY_UNISPG [Options] { bam_list | sample1.bam ...}\n\
"

#define BSIZE 10000 // bundle size

#define MAX_NODE 1000000
#define KMER 31

#define DROP 0.5
#define ERROR_PERC 0.1
#define DBL_ERROR 0.01

#define CHI_WIN 100
#define CHI_THR 50
#define SMALL_EXON 35 // exons smaller than this have a tendency to be missed by long read data

#define IS_FPKM_FLAG 1
#define IS_TPM_FLAG 2
#define IS_COV_FLAG 4

#endif