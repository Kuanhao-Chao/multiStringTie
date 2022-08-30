#include "GArgs.h"
#include "rlink.h"
#include "unispg.h"
#include "helper.h"

void get_fragment_pattern(BundleData* bundle, GList<CReadAln>& readlist, int n, int np, float readcov, GPVec<UnispgGp>** graphs_vec);

void get_read_pattern(int s, float readcov, float rprop, GList<CReadAln>& readlist, int n,GVec<int>& rgno, GVec<int> *rnode, GPVec<UnispgGp>** graphs_vec);