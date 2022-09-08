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
	if(mergeMode) {
	// 	//count_merge_junctions(bundle->readlist,bundle->covflags); // this make sense only if I want to count junctions
		// geneno = build_merge(bundle);

	} 
	// else if (multiMode && (bundle->keepguides.Count() || !eonly)) {
	// 	fprintf(stderr, "This is the multiMiode of the graph.\n");
	// 	count_good_junctions(bundle);
	// 	geneno = build_graphs_multi(bundle, unispg_gp);
	// }
	else if(bundle->keepguides.Count() || !eonly) {
		//fprintf(stderr,"Process %d reads from %lu.\n",bundle->readlist.Count(),bundle->numreads);
		count_good_junctions(bundle);
		geneno = build_graphs_CREATE_UNISPG(bundle, unispg_gp, fidx);
		// geneno = build_graphs_unispg(bundle, unispg_gp);
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