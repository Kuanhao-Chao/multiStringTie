#include "rlink.h"
#include "multist.h"
#include "helper.h"
#include "parse_reads.h"

extern ofstream ratio_file_pos;
extern ofstream ratio_file_neg;

extern bool mergeMode; // "merge" tag. 
extern bool eonly; // "e" tag. for mergeMode includes estimated coverage sum in the merged transcripts
extern int skip_counter;

int infer_transcripts_CREATE_UNISPG(BundleData* bundle, int fidx);
void infer_transcripts_unispg(BundleData* bundle, GPVec<UnispgGp>** graphs_vec);
