#include "rlink.h"
#include "multist.h"
#include "helper.h"
#include "parse_reads.h"

extern ofstream ratio_file_pos;
extern ofstream ratio_file_neg;

extern int skip_counter;

void infer_transcripts_unispg(BundleData* bundle, GPVec<UnispgGp>** graphs_vec);
