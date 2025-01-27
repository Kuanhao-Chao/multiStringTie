#include "infer_tranx_A.h"

void infer_transcripts_APPLY_UNISPG(BundleData* bundle, UnispgGp_APPLY* unispgs) {

	int gene_num = 0;
    int refstart = bundle->start;
    int refend = bundle->end;
	GList<CPrediction>& pred = bundle->pred;

    count_good_junctions_APPLY_UNISPG(bundle);
	
    // geneno = build_graphs_unispg(bundle, unispg);
    // fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
    // fprintf(stderr, "bundle->bpcov[1].Count(): %d\n", bundle->bpcov[1].Count());
    // fprintf(stderr, "bundle->start - end: %d - %d\n", bundle->start, bundle->end);


    /*****************************
     * Iterate through graphs.
     *****************************/
    /*
    { // Debug
        for(int s=0;s<2;s+=1) {
            fprintf(stderr, "\t@@ unispgs[%d]: %d\n", s, unispgs->graph_num[s]);
            for (int i=0; i<unispgs->graph_num[s]; i++) {
                fprintf(stderr, "\t\t unispgs bound: %u - %u \n", unispgs->no2gnode_unispg[s][i][1]->start, unispgs->no2gnode_unispg[s][i][ unispgs->node_nums[s][i]-2 ]->end);
            }
        }
    }
    */

    /*****************************
     ** Variables Definition
     *****************************/
    GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    GPVec<CMTransfrag> *mgt[2]; // merged super-transfrags
    CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag

    for(int s=0;s<2;s++) { // skip neutral bundles -> those shouldn't have junctions
        int g_num = unispgs->graph_num[s];
        transfrag[s] = NULL;
        mgt[s] = NULL;
        tr2no[s] = NULL;

        if (g_num) {
            transfrag[s] = new GPVec<CTransfrag>[g_num];
    		mgt[s] = new GPVec<CMTransfrag>[g_num];
    		GCALLOC(tr2no[s], g_num*sizeof(CTreePat *));

            /*****************************
             * Iterate through graphs to create variable.
             *****************************/
            for (int g=0; g<g_num; g++) {
                // fprintf(stderr, "\t\t (%d) unispgs bound: %u - %u; node_nums: %d\n", s, unispgs->no2gnode_unispg[s][g][1]->start, unispgs->no2gnode_unispg[s][g][ unispgs->node_nums[s][g]-2 ]->end, unispgs->node_nums[s][g]);
                // fprintf(stderr, ">>>>>> s: %d;  g: %d\n", s, g);
                tr2no[s][g]=construct_treepat_APPLY_UNISPG(unispgs->node_nums[s][g], unispgs->gpos[s][g], transfrag[s][g]);

                GStr pattern;
                print_pattern(tr2no[s][g], pattern, unispgs->node_nums[s][g]);
            }
        }
    }

    /*****************************
     **    'get_fragment_pattern_APPLY_UNISPG' function
     **        because of this going throu
     *****************************/
    int global_gidx[2] = {0}; // For tracking what is the last graph that a read overlapped.
    for (int n=0; n < bundle->readlist.Count(); n++) {


        uint read_start = bundle->readlist[n]->start;
        uint read_end = bundle->readlist[n]->end;

        // Only process the reads that are completely in the graph boundary.
        if (read_start >= refstart && read_end <= refend) {


            /*****************************
             * Count fragments in the bundle
             *****************************/
            CReadAln & rd=*(bundle->readlist[n]);

            bundle->frag_len+=rd.len*rd.read_count; 
            double num_frag=rd.read_count;
            for(int i=0;i<rd.pair_idx.Count();i++) {
                // Not counting the fragment if I have seen the pair before
                if(rd.pair_idx[i]>-1 && n>rd.pair_idx[i] && bundle->readlist[rd.pair_idx[i]]->nh) {
                    // only if the read is paired and comes as the first in the pair => count the fragments
                    num_frag-=rd.pair_count[i];
                }
            }
            if(num_frag>epsilon) {
                bundle->num_fragments+=num_frag;
            }

            /*****************************
             * Process the paired-end read.
             *****************************/
            double single_count=rd.read_count;
            for (int i=0; i<bundle->readlist[n]->pair_idx.Count(); i++) {
                // fprintf(stderr, ">> single_count: %d\n", i);
                int np = bundle->readlist[n]->pair_idx[i];
                if (np > -1) {
                    single_count -= bundle->readlist[n]->pair_count[i];
                    // fprintf(stderr, ">> n: %d;   np: %d\n", n, np);
                    if (n < np) {
                        // fprintf(stderr, ">> n < np: %d\n", n < np);
                        get_fragment_pattern_APPLY_UNISPG(bundle, bundle->readlist, n, np, bundle->readlist[n]->pair_count[i], unispgs, global_gidx, transfrag, mgt, tr2no);
                    }
                }
            }
            if (single_count > epsilon) {
                get_fragment_pattern_APPLY_UNISPG(bundle, bundle->readlist, n, -1, single_count, unispgs, global_gidx, transfrag, mgt, tr2no);
            }
        } else if (read_start < refstart && read_end >= refstart) {
            boundary_counter += 1;   
        } else if (read_start <= refend && read_end > refend) {   
            boundary_counter += 1;   
        } else {
            skip_counter += 1;
            skip_counter_nh += bundle->readlist[n]->nh;
            // bundle->readlist[n]->juncs;
            // fprintf(stderr, ">> nh: %d\n", bundle->readlist[n]->nh);
            // bundle->readlist[n]->pair_count;
            // bundle->readlist[n]->pair_idx;
            // fprintf(stderr, ">> read_count: %f\n", bundle->readlist[n]->read_count);
            // bundle->readlist[n]->segs;
            // fprintf(stderr, ">> strand: %d\n", bundle->readlist[n]->strand);

            // fprintf(stderr, "Skipped!!!\n");
        }
    }


    /*****************************
     **    Time to parse the graph!!!
     *****************************/
    for(int s=0;s<2;s++) {
        int g_num = unispgs->graph_num[s];
        for(int g=0;g<g_num;g++) {
            int node_num = unispgs->node_nums[s][g];
            int edge_num = unispgs->edge_nums[s][g];

            /*****************************
             ** 'traverse_dfs' function
             ** 	finished reading bundle -> now create the parents' and children's patterns
             *****************************/
            GVec<bool> visit;
            visit.Resize(node_num);
            GBitVec parents(node_num + edge_num);

            GPVec<CGraphnodeUnispg>** no2gnode = unispgs->no2gnode_unispg;
            GIntHash<int>** gpos = unispgs->gpos;// for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern

            // fprintf(stderr,"traverse graph[%d][%d] now with %d nodes, %d edges and lastgpos=%d....\n",s,g,graphno,edgeno,lastgpos);//edgeno=0;
            traverse_dfs_APPLY_UNISPG(s, g, no2gnode[s][g][0], no2gnode[s][g][node_num-1], parents, node_num, visit, no2gnode, transfrag, edge_num, gpos);
            // fprintf(stderr,"done traversing with edgeno=%d lastgpos=%d\n",edgeno,lastgpos);

            if (node_num) {
                process_transfrags_APPLY_UNISPG(s, node_num, edge_num, no2gnode[s][g], transfrag[s][g],tr2no[s][g],gpos[s][g],pred);

    			GVec<CGuide> guidetrf;
	            GPVec<GffObj>& guides = bundle->keepguides;
	            GVec<int> guidepred; // for eonly keeps the prediction number associated with a guide
    			GVec<int> trflong; // non-redundant long transfrags that I can use to guide the long read assemblies

                gene_num = find_transcripts_APPLY_UNISPG(node_num, edge_num, gpos[s][g], no2gnode[s][g], transfrag[s][g], gene_num, s, guidetrf, guides, guidepred, bundle, trflong);

                // fprintf(stderr, "#######################\n");
                // fprintf(stderr, "## Gene number: %d ##\n", gene_num);
                // fprintf(stderr, "#######################\n");
            }
        }
    }
    // geneno=merge_transfrags(graphno[s][b],no2gnode[s][b], mgt[s][b],gpos[s][b],geneno,s,pred,readlist,guides);

    /*****************************
     ** Output the Stringtie vs redistributed coverage 
     *****************************/   
    // if (graphs_vec[0]->Count() == 0) {
    //     fprintf(stderr, "####### Positive strand!!!\n");
    //     for (int g=0; g<graphs_vec[1]->Count(); g++) {
    //         for (int n=1; n<graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Count()-1; n++) {
    //             float target_cov_ratio = 0.0;
    //             fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[1]->Get(g)->get_refstart(), graphs_vec[1]->Get(g)->get_refend(), g, n, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->nodeid, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end);


    //             float expected_cov_pos = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s.Get(0);
    //             int expected_len_pos = int(graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->len());
    //             float new_cov_pos = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_unispg_s[0];

    //             fprintf(stderr, ">> expected_cov_pos: %f;  new_cov_pos: %f\n", expected_cov_pos/expected_len_pos, new_cov_pos/expected_len_pos);

    //             cov_file_pos << expected_cov_pos << "\t" << new_cov_pos << "\n";
    //             cov_file_pos_norm << (expected_cov_pos/expected_len_pos) << "\t" << (new_cov_pos/expected_len_pos) << "\n";
    //         }
    //     }
    // } else if (graphs_vec[1]->Count() == 0) {
    //     fprintf(stderr, "####### Negative strand!!!\n");
    //     for (int g=0; g<graphs_vec[0]->Count(); g++) {
    //         for (int n=1; n<graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Count()-1; n++) {
    //             float target_cov_ratio = 0.0;
    //             fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[0]->Get(g)->get_refstart(), graphs_vec[0]->Get(g)->get_refend(), g, n, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->nodeid, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end);
                
    //             float expected_cov_neg = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s.Get(0);
    //             int expected_len_neg = int(graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->len());
    //             float new_cov_neg = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_unispg_s[0];

    //             fprintf(stderr, ">> expected_cov_neg: %f;  new_cov_neg: %f.", expected_cov_neg/expected_len_neg, new_cov_neg/expected_len_neg);
                

    //             cov_file_neg << expected_cov_neg << "\t" << new_cov_neg << "\n";
    //             cov_file_neg_norm << (expected_cov_neg/expected_len_neg) << "\t" << (new_cov_neg/expected_len_neg) << "\n";
    //         }
    //     }
    // }
}


void create_graph_param(int s, int g, UnispgGp_APPLY* unispgs, GPVec<CTransfrag> transfrag, GIntHash<int>** gpos, int& lastgpos) {
    // Traverse the graph.
    int node_num = unispgs->no2gnode_unispg[s][g].Count();
    cout << "Nodes in graph (" << s << ", " << g << "): " << unispgs->no2gnode_unispg[s][g].Count() << endl;

    GVec<bool> visit;
	visit.Resize(node_num);
    graph_dfs(s, g, unispgs, unispgs->no2gnode_unispg[s][g][0], visit);
}


void graph_dfs(int s, int g, UnispgGp_APPLY* unispgs, CGraphnodeUnispg* node, GVec<bool>& visit) {
    if (visit[node->nodeid]) {
    } else {
		visit[node->nodeid]=true;
        // fprintf(stderr, "node: %d\n", node->nodeid);
        // fprintf(stderr, "dfs: unispgs[%d][%d] nodeid: %d (%u - %u)\n", s, g, node->nodeid, node->start, node->end);

        for (int c=0; c<node->child.Count(); c++) {
            graph_dfs(s, g, unispgs, unispgs->no2gnode_unispg[s][g][node->child[c]], visit);
        }
    }
}


GBitVec traverse_dfs_APPLY_UNISPG(int s,int g,CGraphnodeUnispg *node,CGraphnodeUnispg *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnodeUnispg> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos){

	// fprintf(stderr,"Traverse node %d;  s: %d;  g: %d;  gno: %d;  edgeno: %d\n",node->nodeid, s, g, gno, edgeno);

	if(visit[node->nodeid]) {
	    // fprintf(stderr,">> node %d visited. \n", node->nodeid);
        // fprintf(stderr, "  Printing parentpat pattern\n");
        // printBitVec(parents);
        // fprintf(stderr, "\n");
		node->parentpat = node->parentpat | parents;
		for(int n=0;n<gno;n++) {
			if(parents[n]) {// add node's children to all parents of node

                // fprintf(stderr, "\tnode->childpat: \n");
                // printBitVec(node->childpat);
                // fprintf(stderr, "\n");
                // fprintf(stderr, "\tno2gnode[s][g][n]->childpat: \n");
                // printBitVec(no2gnode[s][g][n]->childpat);
                // fprintf(stderr, "\n");
				no2gnode[s][g][n]->childpat = no2gnode[s][g][n]->childpat | node->childpat;
            }
			else if(node->childpat[n]) {
                // fprintf(stderr, "\tnode->parentpat: \n");
                // printBitVec(node->parentpat);
                // fprintf(stderr, "\n");
                // fprintf(stderr, "\tno2gnode[s][g][n]->parentpat: \n");
                // printBitVec(no2gnode[s][g][n]->parentpat);
                // fprintf(stderr, "\n");
				no2gnode[s][g][n]->parentpat = no2gnode[s][g][n]->parentpat | node->parentpat;
            }
		}
	}
	else {
		node->childpat.resize(gno+edgeno);
		node->parentpat.resize(gno+edgeno);
		node->parentpat = node->parentpat | parents;
		visit[node->nodeid]=true;
		parents[node->nodeid]=1; // add the node to the parents

    // fprintf(stderr, "  Printing parents pattern\n");
	// printBitVec(parents);
	// fprintf(stderr, "\n");

		if(node->parent.Count()==1 && !node->parent[0]) { // node has source only as parent -> add transfrag from source to node
			GBitVec trpat(gno+edgeno);
			trpat[0]=1;
			trpat[node->nodeid]=1;

			int key=edge(0,node->nodeid,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			GVec<int> nodes;
			nodes.cAdd(0);
			nodes.Add(node->nodeid);
			CTransfrag *tr=new CTransfrag(nodes,trpat,trthr);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add source transfrag[%d][%d]= %d and pattern",s,g,transfrag[s][g].Count());
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
            */

			transfrag[s][g].Add(tr);
		}

		int n=node->child.Count();
		if(node != sink && !n) {
			node->child.Add(sink->nodeid);  // add sink to the node's children
			sink->parent.Add(node->nodeid); // add node to sink's parents
			// create the transfrag that ends the node
			GBitVec trpat(gno+edgeno);
			trpat[node->nodeid]=1;
			trpat[gno-1]=1;

			int key=edge(node->nodeid,gno-1,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;

			GVec<int> nodes;
			nodes.Add(node->nodeid);
			nodes.Add(sink->nodeid);
			CTransfrag *tr=new CTransfrag(nodes,trpat,trthr);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add sink transfrag[%d][%d]= %d for nodeid=%d and pattern:",s,g,transfrag[s][g].Count(),node->nodeid);
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
            */

			transfrag[s][g].Add(tr);
			n++;
	    }
		/*
		fprintf(stderr,"Add %d children of node %d (%d-%d): ",n,node->nodeid,node->start,node->end);
		for(int i=0;i<n;i++) fprintf(stderr," %d",node->child[i]);
		fprintf(stderr,"\n");
		*/

		//edgeno+=n; // this will have to be deleted in the end; now I put it so that I can check equivalence with the one computed when creating the graph

	    for(int i=0; i< n; i++) { // for all children
            // fprintf(stderr, "  i: %d\n", i);
	    	GBitVec childparents=parents;

    // fprintf(stderr, "  Printing childparents pattern\n");
	// printBitVec(childparents);
	// fprintf(stderr, "\n");

	    	int min=node->nodeid; // nodeid is always smaller than child node ?
	    	int max=node->child[i];
	    	if(min>max) {
	    		max=node->nodeid; // nodeid is always smaller than child node ?
	    		min=node->child[i];
	    	}

    // fprintf(stderr, "  Printing min-max: %d-%d\n", min, max);

			int key=edge(min,max,gno);
			int *pos=gpos[s][g][key];

			if(pos!=NULL) {
				childparents[*pos]=1; // add edge from node to child to the set of parents from child
				node->childpat[*pos]=1; // add edge from node to child to the set of node children
			}
			// fprintf(stderr,"Call for child %d with id=%d\n",node->child[i],no2gnode[s][g][node->child[i]]->nodeid);
	    	node->childpat = node->childpat | traverse_dfs_APPLY_UNISPG(s,g,no2gnode[s][g][node->child[i]],sink,childparents,gno,visit,no2gnode,transfrag,edgeno,gpos);
	    }
	} // end else from if(visit[node->nodeid])

	GBitVec children = node->childpat;
    // fprintf(stderr, "  Printing children pattern\n");
	// printBitVec(children);
	// fprintf(stderr, "\n");
	children[node->nodeid]=1;

	return(children);
}
