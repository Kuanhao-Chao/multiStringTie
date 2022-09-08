#ifndef __UNISPG_A_H__
#define __UNISPG_A_H__
#pragma once

#include "global_params_A.h"
#include "../unispg.h"

#include <unordered_map>

struct UnispgGp_APPLY:public UnispgGp {
	UnispgGp_APPLY(int refstart_i, int refend_i) {
		refstart = refstart_i;
		refend = refend_i;
		for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
			int s=sno/2; // adjusted strasnd due to ignoring neutral strand
			no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];

			// transfrag_unispg[s] = new GPVec<CTransfrag>[20000];
			// source_gp[s] = new CGraphnodeUnispg[1];
			// sink_gp[s] = new CGraphnodeUnispg[1];
		}
	}
};
#endif

