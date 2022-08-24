#include "GArgs.h"
#include "rlink.h"
#include "unispg.h"

void get_fragment_pattern(GList<CReadAln>& readlist, int n, int np, GPVec<UnispgGp>** graphs_vec);

void get_read_pattern(int s, GList<CReadAln>& readlist, int n,GVec<int> rgno, GVec<int> *rnode, GPVec<UnispgGp>** graphs_vec);