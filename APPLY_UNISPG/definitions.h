
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