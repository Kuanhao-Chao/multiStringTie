#include "infer_tranx_C.h"
int infer_transcripts_CREATE_UNISPG(BundleData* bundle, UnispgGp_CREATE* unispg_gp, int fidx) {
    int geneno=0;
	//DEBUG ONLY: 	showReads(refname, readlist);

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
	if(bundle->keepguides.Count() || !eonly) {
		//fprintf(stderr,"Process %d reads from %lu.\n",bundle->readlist.Count(),bundle->numreads);
		count_good_junctions(bundle);
		geneno = build_graphs_CREATE_UNISPG(bundle, unispg_gp, fidx);
		// geneno = build_graphs_unispg(bundle, unispg_gp);

		// fprintf(stderr, "#######################\n");
		// fprintf(stderr, "## Gene number: %d ##\n", geneno);
		// fprintf(stderr, "#######################\n");
	}

/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(e):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
	return(geneno);
}
