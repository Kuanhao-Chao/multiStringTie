#include "parse_reads_A.h"

void get_fragment_pattern_APPLY_UNISPG(BundleData* bundle, GList<CReadAln>& readlist, int n, int np, float readcov, UnispgGp_APPLY* unispgs, int* global_gidx, GPVec<CTransfrag> **transfrag, GPVec<CMTransfrag>** mgt, CTreePat ***tr2no) {
    // fprintf(stderr, ">> get_fragment_pattern_APPLY_UNISPG is called. \n");

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
	 * Step 1: Calculating the redistribution ratio (for each read).
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

        /*****************************
         * Step 1-1-1: main read!
         *****************************/
    	// rgno: It stores the bundle<->graph indices in a BundleData that a read belongs to.
		int rgno = 0;
	    // rnode: It stores the node indices in the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> rnode;
        // Read pattern
        GBitVec rpat;
		if(readlist[n]->nh) {
			get_read_pattern_APPLY_UNISPG(s, readcov, rprop[s], readlist, n, rnode, unispgs, rpat, global_gidx, transfrag, mgt, tr2no);
        }

        /*****************************
         * Step 1-1-2: paired read!
         *****************************/
    	// pgno: It stores the paired bundle<->graph indices in a BundleData that a read belongs to.
		int pgno = 0;
	    // pnode: It stores the paired node indices in the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> pnode;
        // Paired pattern
        GBitVec ppat;
        if(np>-1 && readlist[np]->nh) {
			get_read_pattern_APPLY_UNISPG(s, readcov, rprop[s], readlist, np, pnode, unispgs, ppat, global_gidx, transfrag, mgt, tr2no);
		}

        int usedp=0;

        // set read pattern
        rpat.clear();
        int i=0;
        int clear_rnode=0;

    }
}


/*****************************
 * (1) confirm each read belongs to which nodes
 * (2) the node coverage gonna be updated
 *****************************/
void get_read_pattern_APPLY_UNISPG(int s, float readcov, float rprop, GList<CReadAln>& readlist, int n, GVec<int>& nodes, UnispgGp_APPLY* unispgs, GBitVec& pat, int* global_gidx, GPVec<CTransfrag> **transfrag, GPVec<CMTransfrag>** mgt, CTreePat ***tr2no) {
	// nodes: It stores the node indices in the bundle<->graph indices in a BundleData that a read belongs to.
	// fprintf(stderr, ">> Inside 'get_read_pattern_APPLY_UNISPG'\n");
	int lastgnode=-1;
	int lastngraph=-1;

	/*****************************
	 * Definition: Read information.
	 *****************************/
	int ncoord=readlist[n]->segs.Count();
	uint r_start = readlist[n]->start;
	uint r_end = readlist[n]->end;
	// int k=0; // need to keep track of coordinates already added to coverages of graphnodes
	int kmer=KMER-1; //f1
	bool gn_r_overlap = false;

	GIntHash<bool> hashnode;
    // Here, I need to link a read to the node/nodes that it belongs to
	/*****************************
	 * Using the whole read start/end to anchor the graph.
	 * (1) Iterate through graphs
	 * (2) Check if a read overlaps the graph.
	 * (3) Iterate through segments in a read & nodes in the graph
	 *     (3.1) Check if a read segment overlap a node.
	 *     (3.2) Create a transfrag for a read.
	 *****************************/

	// fprintf(stderr, ">> Read (%d - %d)\n", r_start, r_end);

	int g_sidx = 0;
	if (global_gidx[s] == 0) {
		g_sidx = 0;
	} else {
		g_sidx = global_gidx[s]-1;
	}
	
	/*****************************
	 * (1) Iterate through graphs
	 *****************************/
	int g = 0; 
	int node_num = 0;
	int edge_num = 0;
	// break when gn_r_overlap equals true.
	for (g=g_sidx; (g<unispgs->graph_num[s] && !gn_r_overlap); g++) {
		// We can create different transfrag based on different graphs.
		/*****************************
		 * Definition: Graph information.
		 *****************************/
		// fprintf(stderr, ">> g_local: %d\n", g);
		int g_start = unispgs->no2gnode_unispg[s][g][1]->start;
		int g_end = unispgs->no2gnode_unispg[s][g][ unispgs->node_nums[s][g]-2 ]->end;
		int g_node_num = unispgs->node_nums[s][g];

		/*****************************
	 	 * (2) Check if a read overlaps the graph.
		 *****************************/
        bool g_r_overlap = segs_overlap(g_start, g_end, r_start, r_end);
		if (g_r_overlap) {
			/*****************************
			 * (2.1) We need to create a pattern for this read based on the overlapped graph.
			 *****************************/
			pat.reset();

			node_num = unispgs->node_nums[s][g];
			edge_num = unispgs->edge_nums[s][g];
			// fprintf(stderr, ">> node_num: %d;  edge_num: %d\n", node_num, edge_num);
			pat.resize(node_num + edge_num);
			/*****************************
	 		 * (3) Iterate through segments in a read & nodes in the graph
			 *****************************/
			int nidx = 1; 
			int k = 0;
			while (nidx < g_node_num-1 && k < ncoord) {
				// fprintf(stderr, ">> s: %d; n: %d; global_gidx[0]: %d; global_gidx[1]: %d; k: %d\n", s, n, global_gidx[0], global_gidx[1], k);

				/*****************************
				 * (3.1) rseq-related information.
				 *****************************/
				uint rseg_start = readlist[n]->segs[k].start;
				uint rseg_end = readlist[n]->segs[k].end;
				// fprintf(stderr, ">> rseg (%u - %u)\n", rseg_start, rseg_end);

				/*****************************
				 * (3.2) node-related information
				 *****************************/
				CGraphnodeUnispg *node = unispgs->no2gnode_unispg[s][g][nidx];
				uint n_start = node->start;
				uint n_end = node->end;
				// fprintf(stderr, ">> node (%u - %u)\n", n_start, n_end);

				/*****************************
				 * (3.3) Check if rseq * node are overlapped.
				 *****************************/

    			bool intersect=false;
				while (k < ncoord) {
					// fprintf(stderr, ">> k: %d,  ncoord: %d\n", k, ncoord);
					rseg_start = readlist[n]->segs[k].start;
					rseg_end = readlist[n]->segs[k].end;
					// fprintf(stderr, ">> rseg (%u - %u)\n", rseg_start, rseg_end);
					int bp = node->calOverlapLen(rseg_start, rseg_end);
					if(bp) {
						intersect = true;
						gn_r_overlap = true;
						// fprintf(stderr, ">> bp: %d \n", bp);
						unispgs->no2gnode_unispg[s][g][nidx]->cov_unispg_s[0] += rprop*bp*readcov;
						// fprintf(stderr, ">> node->cov_s[0]: %f\n", unispgs->no2gnode_unispg[s][g][nidx]->cov_unispg_s[0]);
						if (rseg_end <= n_end) k++;
						else break; 
					} else break;
				}

				if (intersect) {
					/*****************************
					 ** (3.4) Overlap occurs => set the edge bit first. 
					 *****************************/
					int *pos;

					if (lastgnode == -1) {
						// fprintf(stderr, "nodes.Add(nidx): nidx: %d\n", nidx);
						nodes.Add(nidx);
						/*****************************
						 ** (3.5) Overlap occurs => set the node bit 
						*****************************/
						pat[nidx] = 1;
						lastgnode = nidx;
					} else if (lastgnode > -1) {
						int min = lastgnode;
						int max = nidx; 	
						pos = unispgs->gpos[s][g][edge(min,max,g_node_num)];
						// fprintf(stderr, "### @@ Adding edge (%d - %d);  g_node_num: %d\n", min, max, g_node_num);
						// fprintf(stderr, ">>> pos: %d\n", *pos);
						if(pos!=NULL) {
							// fprintf(stderr, "nodes.Add(nidx): nidx: %d\n", nidx);
							nodes.Add(nidx);
							/*****************************
							 ** (3.5) Overlap occurs => set the node bit 
							*****************************/
							pat[nidx] = 1;
							pat[*pos] = 1;
						// 	// fprintf(stderr, "### Adding edge (%d - %d): %d \n", min, max, *pos);
							lastgnode = nidx;
						} else {
							//  I am discarding this node!!!
						}
					}
					/*****************************
					 ** (3.4.1) First intuition!! Skip adding this node if the edge does not exist
					 *****************************/
					// if(pos!=NULL) {
					// 	fprintf(stderr, "nodes.Add(nidx): nidx: %d\n", nidx);
					// 	nodes.Add(nidx);
					// 	/*****************************
					// 	 ** (3.5) Overlap occurs => set the node bit 
					// 	*****************************/
					// 	pat[nidx] = 1;
					// 	lastgnode = nidx;
					// }
				}
				nidx++;
			}

		} else if (g_end <= r_start) {
			// Move on to the next graph!!
			// fprintf(stderr, "|(g)-------(g)|   (r).......(r)\n");
			global_gidx[s] += 1;
			// fprintf(stderr, ">> s: %d; global_gidx[s]: %d; g_local: %d\n", s, global_gidx[s], g);
			continue;
		} else {
			// (r).......(r)  |(g)-------(g)|
			// It's done processing this read. => break!
			break;
		}
	}

	g -= 1;
	if (gn_r_overlap) {
		// fprintf(stderr, ">> Printing bitvector: ");
		// printBitVec(pat);
		// fprintf(stderr, "\n");

      	/*****************************
         * (4) Update the abundance of transfrag 
		 *   (4.1) Update the abundance.
		 *   (4.2) push the read pattern into transfrag
		 *   (4.3) Iteratively create the CTreePat.
         *****************************/
		CTransfrag *t=update_read_pattern_abund_APPLY_UNISPG(s, g, node_num, unispgs->gpos[s][g], pat, readcov, nodes, transfrag, mgt, tr2no);


      	// /*****************************
        //  * (4) Update the MTransfrag
        //  *****************************/
		// int mgt_idx = 0;
		// while (mgt_idx < mgt[s][g].Count()) {
		// 	if(mgt[s][g][mgt_idx]->transfrag == t) {
		// 		for(int x=0;x<readlist[n]->pair_idx.Count();x++) {
		// 			mgt[s][g][mgt_idx]->read.Add(readlist[n]->pair_idx[x]);
		// 		}
		// 		//mgt[s][rgno[s]][i]->read.Add(n);
		// 		mgt[s][g][mgt_idx]->len=readlist[n]->len;
		// 		break;
		// 	}
		// 	mgt_idx++;
		// }
		// if(mgt_idx == mgt[s][g].Count()) { // MTransfrag not found
		// 	CMTransfrag *mt=new CMTransfrag(t);
		// 	for(int x=0;x<readlist[n]->pair_idx.Count();x++) {
		// 		mt->read.Add(readlist[n]->pair_idx[x]);
		// 	}
		// 	mt->len=readlist[n]->len;
		// 	mgt[s][g].Add(mt);
		// }
	}
}


CTransfrag* update_read_pattern_abund_APPLY_UNISPG(int s,int g, int node_num, GIntHash<int>& gpos, GBitVec& pat, float abundance, GVec<int>& nodes, GPVec<CTransfrag> **transfrag, GPVec<CMTransfrag>** mgt, CTreePat ***tr2no){

	// fprintf(stderr, ">> Inside 'update_read_pattern_abund_APPLY_UNISPG'!!\n");
	/*
	{ // DEBUG ONLY
		// fprintf(stderr,"Update transfrag[%d][%d] longread=%d:",s,g,is_lr);
		for(int i=0;i<nodes.Count();i++) fprintf(stderr," %d",nodes[i]);
		fprintf(stderr," with abundance=%f in graph[%d][%d]\n",abundance,s,g);
	}
	*/
	CTransfrag *t = findtrf_in_treepat_APPLY_UNISPG(node_num, gpos, nodes,pat,tr2no[s][g]);
	if(!t) { // t is NULL
	  t=new CTransfrag(nodes,pat,0);

	  /*
		{ // DEBUG ONLY
			fprintf(stderr,"Add update transfrag[%d][%d]=%d and pattern",s,g,transfrag[s][g].Count());
			fprintf(stderr," (check nodes:");
			for(int i=0;i<nodes.Count();i++) fprintf(stderr," %d",nodes[i]);
			fprintf(stderr,")");
			//printBitVec(pattern);
			fprintf(stderr,"\n");
		}
	   */

	  transfrag[s][g].Add(t);

	  // nodes.Sort() : nodes should be sorted; if they are not then I should update to sort here
	  CTreePat *tree=tr2no[s][g];
	  for(int n=0;n<nodes.Count();n++) {
		  CTreePat *child;
		  if(n) { // not the first node in pattern
			  int *pos=gpos[edge(nodes[n-1],nodes[n],node_num)];
			  if(pos && pat[*pos]) // there is an edge between nodes[n-1] and nodes[n]
				  child=tree->settree(node_num-1-nodes[n-1]+nodes[n]-nodes[n-1]-1,nodes[n],2*(node_num-nodes[n]-1));
			  else child=tree->settree(nodes[n]-nodes[n-1]-1,nodes[n],2*(node_num-nodes[n]-1));
		  }
		  else child=tree->settree(nodes[n]-1,nodes[n],2*(node_num-nodes[n]-1));
		  tree=child;
	  }
	  tree->tr=t;
	}

	//fprintf(stderr,"Set short read\n");
	t->shortread=true;
	t->abundance+=abundance;

	return(t);
}