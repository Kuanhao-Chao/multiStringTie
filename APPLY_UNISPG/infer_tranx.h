#include "rlink.h"
#include "multist.h"
#include "helper.h"
#include "parse_reads.h"

extern ofstream ratio_file;

void infer_transcripts_unispg(BundleData* bundle, GPVec<UnispgGp>** graphs_vec);
