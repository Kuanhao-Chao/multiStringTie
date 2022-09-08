#include "parse_reads_A.h"

void get_fragment_pattern_APPLY_UNISPG(BundleData* bundle, GList<CReadAln>& readlist, int n, int np, float readcov, GPVec<UnispgGp_APPLY>** graphs_vec, int* global_gidx) {
    fprintf(stderr, "get_fragment_pattern_APPLY_UNISPG is called. \n");

	int refstart = bundle->start;
    int refend = bundle->end;
	uint rstart=readlist[n]->start;
	uint rend=readlist[n]->end;
	if(np>-1 && readlist[np]->end>rend) {
		rend=readlist[np]->end;
	}
	
	// The neg / pos proportion
	float rprop[2]={0,0}; // by default read does not belong to any strand

	/*****************************
	 * Step 1: Calculatingt the redistribution ratio (for each read).
	 *****************************/
	// nh => number of reported alignments that contain the query in the current record.
	if(readlist[n]->nh && !readlist[n]->strand && np>-1 && readlist[np]->nh && !readlist[np]->strand) {
		/*****************************
		 *  both reads are unstranded
		 *****************************/
		redistribute_unstranded_rcov(rprop, bundle->bpcov, refstart, refend, rstart, rend);
		// fprintf(stderr, "rprop[0]: %f; rprop[1]: %f\n", rprop[0], rprop[1]);
	} else {
		/*****************************
		 *  Any of the read is stranded.
		 *****************************/
		if(readlist[n]->nh) {
			/*****************************
			 *  The number of reported alignments of current read > 0
			 *****************************/
			if(!readlist[n]->strand) { // the paired read is not present otherwise it would have the same strand
				redistribute_unstranded_rcov(rprop, bundle->bpcov, refstart, refend, rstart, rend);
				// fprintf(stderr, "rprop[0]: %f; rprop[1]: %f\n", rprop[0], rprop[1]);
			} else {
				if(readlist[n]->strand==-1) rprop[0]=1;
				else rprop[1]=1;
			}
		} else if(np>-1 && readlist[np]->nh) {
			/*****************************
			 *  The number of reported alignments of the pair of the current read > 0
			 *****************************/
			if(!readlist[np]->strand) { // the paired read is not present otherwise it would have the same strand
				redistribute_unstranded_rcov(rprop, bundle->bpcov,refstart, refend, rstart, rend);
				// fprintf(stderr, "rprop[0]: %f; rprop[1]: %f\n", rprop[0], rprop[1]);
			} else {
				if(readlist[np]->strand==-1) rprop[0]=1;
				else rprop[1]=1;
			}
		}
	}

	/*****************************
	 * Step 1: Generate the read pattern !!!
	 *  both reads are unstranded
	 *****************************/
	for(int s=0;s<2;s++) {

        /*****************************
         * Step 1-1: After the `get_read_pattern_APPLY_UNISPG`, 
         *   (1) confirm each read belongs to which nodes
         *   (2) the node coverage gonna be updated
         *****************************/
    	// rgno: It stores the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> rgno;
	    // rnode: It stores the node indices in the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> *rnode=new GVec<int>[readlist[n]->segs.Count()];

		if(readlist[n]->nh) {
			get_read_pattern_APPLY_UNISPG(s, readcov, rprop[s], readlist, n, rgno, rnode, graphs_vec, global_gidx);
        }

    	// pgno: It stores the paired bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> pgno;
	    // pnode: It stores the paired node indices in the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> *pnode=NULL;

        if(np>-1 && readlist[np]->nh) {
			pnode=new GVec<int>[readlist[np]->segs.Count()];
			get_read_pattern_APPLY_UNISPG(s, readcov, rprop[s], readlist, np, rgno, rnode, graphs_vec, global_gidx);
		}


        /*****************************
         * Step 1-2: actually create a pattern for a read
         *****************************/
        // Read pattern
        GBitVec rpat;
        // Paired pattern
        GBitVec ppat;
        int usedp=0;

        // set read pattern
        rpat.clear();
        // rpat.resize(graphno[s][rgno[r]]+edgeno[s][rgno[r]]);
        int i=0;
        int clear_rnode=0;

    }
}


/*****************************
 * (1) confirm each read belongs to which nodes
 * (2) the node coverage gonna be updated
 *****************************/
void get_read_pattern_APPLY_UNISPG(int s, float readcov, float rprop, GList<CReadAln>& readlist, int n,GVec<int> &rgno, GVec<int> *rnode, GPVec<UnispgGp_APPLY>** graphs_vec, int* global_gidx) {

	// rgno: It stores the bundle<->graph indices in a BundleData that a read belongs to.
	// rnode: It stores the node indices in the bundle<->graph indices in a BundleData that a read belongs to.

	fprintf(stderr, ">> Inside 'get_read_pattern_APPLY_UNISPG'\n");
	int lastgnode=-1;
	int lastngraph=-1;
	int ncoord=readlist[n]->segs.Count();
	// int k=0; // need to keep track of coordinates already added to coverages of graphnodes
	int kmer=KMER-1; //f1

	GIntHash<bool> hashnode;

    // Here, I need to link a read to the node/nodes that it belongs to

	/*****************************
	 * Using the whole read start/end to anchor the graph.
	 * (1) Iterate through segments in nodes. Check which node it belongs to. 
	 * (2) Iterate through graphs
	 * (3) Iterate through nodes in the graph
	 *****************************/
	uint r_start = readlist[n]->start;
	uint r_end = readlist[n]->end;
	fprintf(stderr, ">> Read (%d - %d)\n", r_start, r_end);

	for (int k=0; k<ncoord; k++) {
    	fprintf(stderr, ">> s: %d; n: %d; global_gidx[0]: %d; global_gidx[1]: %d; k: %d\n", s, n, global_gidx[0], global_gidx[1], k);


		uint rseg_start = readlist[n]->segs[k].start;
		uint rseg_end = readlist[n]->segs[k].end;
		for (int g=global_gidx[s]; g<graphs_vec[s]->Count(); g++) {
		// for (int g=0; g<graphs_vec[s]->Count(); g++) {
			fprintf(stderr, ">> g_local: %d\n", g);
			int g_start = graphs_vec[s]->Get(g)->get_refstart();
			int g_end = graphs_vec[s]->Get(g)->get_refend();


			if (g_end < r_start) {
				fprintf(stderr, "|(g)-------(g)|   (r).......(r)\n");
				global_gidx[s] += 1;
				fprintf(stderr, ">> s: %d; global_gidx[s]: %d; g_local: %d\n", s, global_gidx[s], g);
				continue;
			}

            bool g_rseq_overlap = segs_overlap(g_start, g_end, rseg_start, rseg_end);

			if (g_rseq_overlap) {


				// Iterating through nodes in the graph
				for (int nidx=1; nidx<graphs_vec[s]->Get(g)->no2gnode_unispg[s][0].Count()-1; nidx++) {
					// fprintf(stderr, ">> graphs_vec[s]->Get(g)->no2gnode_unispg[s][0].Count(): %d\n", graphs_vec[s]->Get(g)->no2gnode_unispg[s][0].Count());

					CGraphnodeUnispg *node = graphs_vec[s]->Get(g)->no2gnode_unispg[s][0][nidx];
					// fprintf(stderr, ">> node: %p\n", node);
					// fprintf(stderr, ">> node: %u - %u\n", node->start, node->end);

					uint n_start = node->start;
					uint n_end = node->end;
					// fprintf(stderr, ">> graphs_vec[s]->Get(g)->no2gnode_unispg[s][0]: %u - %u\n", n_start, n_end);

					// int bp = readlist[n]->segs[k].overlapLen(n_start, n_end);
					int bp = node->calOverlapLen(rseg_start, rseg_end);

					// int bp = calOverlapLen(n_start, n_end, rseg_start, rseg_end);
					// // // int bp = readlist[n]->segs[k].overlapLen();
					// // // int bp = overlapLen(n_start, n_end, rseg_start, rseg_end);
					if(bp) {
						// intersect=true;
						fprintf(stderr, ">> bp: %d \n", bp);
						/*****************************
						 ** Update the realist if it is not a unitig.
						 *****************************/
						// node->cov_s[0]+=rprop*bp*readcov;
						
						graphs_vec[s]->Get(g)->no2gnode_unispg[s][0][nidx]->cov_unispg_s[0] += rprop*bp*readcov;

						fprintf(stderr, ">> node->cov_s[0]: %f\n", graphs_vec[s]->Get(g)->no2gnode_unispg[s][0][nidx]->cov_unispg_s[0]);

						// fprintf(stderr, ">> node->cov_s[0]: %f; (rprop: %f; bp: %d; readcov: %f) \n", rprop*bp*readcov, rprop, bp, readcov);


						// if(readlist[n]->segs[k].end<=n_end) ;
						// else break;

						// else if (n_end < rseg_start) {
						// 	// |(n)-------(n)|  (rseg).......(rseg)
						// 	continue;
						// }
					}
					if (rseg_end <= n_end) {
						// (rseg).......(rseg)  |(n)-------(n)|
						break;
					} 
				}

			} else if (rseg_end < g_start) {
				// (rseg).......(rseg)  |(g)-------(g)|
				fprintf(stderr, "(rseg).......(rseg)  |(g)-------(g)|\n");
				break;
			} else if (g_end < rseg_start) {
				// |(g)-------(g)|  (rseg).......(rseg)
				fprintf(stderr, "|(g)-------(g)|  (rseg).......(rseg)\n");
				continue;
			}
		}
	}
}