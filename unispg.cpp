#include "unispg.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

extern UnispgGp* unispg_gp;
extern GVec<int> current_gidx;
extern FILE* uinigraph_out;

void UnispgGp::ProcessSample(GStr sample_name) {
    samples.Add(sample_name);
    for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
        int s=sno/2; // adjusted strand due to ignoring neutral strand
        current_gidx[s] = 0;
    }
}

void UnispgGp::AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode, int lclg_limit) {
    // Next lclg: no2gnode+1 !!    


    int sample_num = unispg_gp->samples.Count();
    if (fidx == 0) {
        for(int g=0; g<lclg_limit; g++) {
            if(no2gnode[g].Count()) {
                fprintf(stderr, "\n*****************************\n");
                fprintf(stderr, "*********** AddGraph ********\n");
                fprintf(stderr, "*****************************\n");
                int unispg_idx = current_gidx[s];
                // This is the first unispg. Simply add copy it into UnispgGp
                for (int i=0; i<no2gnode[g].Count(); i++) {
                    CGraphnode* node = new CGraphnode(no2gnode[g][i]);

                    GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                    GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                    GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);

                    CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(sample_num, node->start, node->end, node->nodeid, is_passed_s, cov_s, capacity_s, true, node->cov, node->capacity, node->rate);

                    // no2gnodeGp_unispg[s][unispg_idx].Add(node);
                    no2gnode_unispg[s][unispg_idx].Add(node_unispg);
                }

                // Linking parent and child
                for (int i=0; i<no2gnode[g].Count(); i++) {
                    fprintf(stderr, "CGraphnode parent: ");
                    for(int p=0; p<no2gnode[g][i]->parent.Count(); p++) {
                        fprintf(stderr, "%d  ", no2gnode[g][i]->parent[p]);
                        fprintf(stderr, "%d  ", no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->parent[p] ]->nodeid);

                        no2gnode_unispg[s][unispg_idx][i]->parent.Add(no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->parent[p] ]);
                    }
                    fprintf(stderr, "\nCGraphnode child: ");
                    for(int c=0; c<no2gnode[g][i]->child.Count(); c++) {
                        fprintf(stderr, "%d  ", no2gnode[g][i]->child[c]); 
                        fprintf(stderr, "%d  ", no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->child[c] ]->nodeid); 

                        no2gnode_unispg[s][unispg_idx][i]->child.Add(no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->child[c] ]);
                    }
                    fprintf(stderr, "\n");
                }





                GStr node_g(g);
                GStr strand_symbol;
                if (s == 0) {
                    strand_symbol = "-";
                } else if (s == 1) {
                    strand_symbol = "+";
                }
                GStr bundle_start("");
                GStr bundle_end("");

                /****************
                 **  Writing out the visualization graph for the local graph.
                ****************/
                bundle_start = int(no2gnode[g][1]->start);
                bundle_end = int(no2gnode[g][no2gnode[g].Count()-2]->end);
                fprintf(stderr,"Traversing the created graph!!!\n");
                // fprintf(stderr,"Digraph %d_%d_%d_%d {", bdata->start,bdata->end, s, g);
                for(int nd=0;nd<no2gnode[g].Count();nd++) {
                    fprintf(stderr,"%d[start=%d end=%d cov=%f];",nd,no2gnode[g][nd]->start,no2gnode[g][nd]->end,no2gnode[g][nd]->cov);
                    // exon_tmp.clear();
                    // exon_tmp.push_back(no2gnode[g][nd]->start);
                    // exon_tmp.push_back(no2gnode[g][nd]->end);
                    // // exon_tmp.push_back(nd*3);
                    // // exon_tmp.push_back(nd*3+1);
                    // exonIntervals.push_back(exon_tmp);
                    // for (int i = no2gnode[g][nd]->start; i < no2gnode[g][nd]->end; i++) {
                    // 	fprintf(node_cov_bed, "chr22\t%d\t%d\tNODE\t%f\t%s\n", i, i+1, no2gnode[g][nd]->cov, strand_symbol.chars());
                    // }
                    int node_start = 0;
                    int node_end = 0;
                    GStr node_nd(nd);
                    GStr node_name = "Node_" + bundle_start + "_" + bundle_end + "_ " + node_g + "_" + node_nd;
                    fprintf(stderr, "node_name: %s\n", node_name.chars());

                    if (nd == 0) {
                        node_start = no2gnode[g][1]->start-50;
                        node_end = no2gnode[g][1]->start;
                    } else if (nd == no2gnode[g].Count()-1){
                        node_start = no2gnode[g][no2gnode[g].Count()-2]->end-1;
                        node_end = 	no2gnode[g][no2gnode[g].Count()-2]->end+50;
                    } else {
                        node_start = no2gnode[g][nd]->start-1;
                        node_end = no2gnode[g][nd]->end;		
                    }


                    if (nd == 0) {
                        if(s == 0) {
                            fprintf(node_cov_neg_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%f\t+\n", node_start, node_end, node_name.chars(), 0);
                        } else if (s == 1) {
                            fprintf(node_cov_pos_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%f\t-\n", node_start, node_end, node_name.chars(), 0);
                        }
                    } else if (nd == no2gnode[g].Count()-1){
                        if(s == 0) {
                            // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\tNODE\t%f\t+\n", no2gnode[g][no2gnode[g].Count()-2]->end, no2gnode[g][no2gnode[g].Count()-2]->end+200, 0);
                            fprintf(node_cov_neg_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%f\t+\n", node_start, node_end, node_name.chars(), 0);
                        } else if (s == 1) {
                            // fprintf(node_cov_pos_bed, "chr22\t%d\t%d\tNODE\t%f\t-\n", no2gnode[g][no2gnode[g].Count()-2]->end, no2gnode[g][no2gnode[g].Count()-2]->end+200, 0);

                            fprintf(node_cov_pos_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%f\t-\n", node_start, node_end, node_name.chars(), 0);
                        }
                    } else {
                        if(s == 0) {
                            // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\t%f\t%s\n", no2gnode[g][nd]->start, no2gnode[g][nd]->end, no2gnode[g][nd]->cov, strand_symbol.chars());
                            fprintf(node_cov_neg_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%f\t%s\n", node_start, node_end, node_name.chars(), no2gnode[g][nd]->cov, strand_symbol.chars());
                        } else if (s == 1) {
                            fprintf(node_cov_pos_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%f\t%s\n", node_start, node_end, node_name.chars(), no2gnode[g][nd]->cov, strand_symbol.chars());
                        }
                    }
                }
                for(int nd=0;nd<no2gnode[g].Count()-1;nd++) {
                    // fprintf(stderr,"Node %d with parents:",i);
                    GStr node_parent_nd(nd);
                    for(int c=0;c<no2gnode[g][nd]->child.Count();c++) {
                        GStr node_child_nd(no2gnode[g][nd]->child[c]);
                        fprintf(stderr,"%d->",nd);			
                        fprintf(stderr,"%d;",no2gnode[g][nd]->child[c]);
                        GStr junction_name = "Junc_" + bundle_start + "_" + bundle_end + "_" + node_g + "_" + node_parent_nd + "->" + node_child_nd;
                        fprintf(stderr, "junction_name: %s\n", junction_name.chars());
                        
                        int junc_start = 0;
                        int junc_end = 0;
                        if (nd == 0) {
                            // It's the source node.
                            junc_start = no2gnode[g][1]->start;
                        } else {
                            junc_start = no2gnode[g][nd]->end;
                        }
                        if (no2gnode[g][ no2gnode[g][nd]->child[c] ] -> start == 0) {
                            // The node goes to the sink.
                            junc_end = no2gnode[g][no2gnode[g].Count()-2]->end;
                        } else {
                            junc_end = no2gnode[g][ no2gnode[g][nd]->child[c] ] -> start;
                        }
                        if(s == 0) {
                            fprintf(edge_cov_neg_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%d\t%s\n", junc_start, junc_end, junction_name.chars(), 10, strand_symbol.chars());
                        } else if (s == 1) {
                            fprintf(edge_cov_pos_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%d\t%s\n", junc_start, junc_end, junction_name.chars(), 10, strand_symbol.chars());
                        }
                    }
                }
                fprintf(stderr,"}\n");



                /****************
                 **  Writing out the visualization graph for the global graph.
                ****************/
                fprintf(stderr, "&& unispg_gp->current_gidx: %d\n", unispg_gp->current_gidx[s]-1);
                // GPVec<CGraphnode>** unispg_gp->no2gnode_unispg = unispg_gp->get_no2gnodeGp();
                fprintf(stderr, "unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count(): %d \n", unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count());
                // int refstart_unispg = unispg_gp->no2gnode_unispg[s][g][1]->start;
                // int refend_unispg = unispg_gp->no2gnode_unispg[s][g][unispg_gp->no2gnode_unispg[s][g].Count()-2]->end;
                // GVec<int>* graphno_unispg = unispg_gp->get_graphnoGp();
                // GVec<int>* edgeno_unispg = unispg_gp->get_edgenoGp();

                fprintf(stderr,"Traversing the universal splice graph!!!\n");
                // fprintf(stderr,"Unispg %d_%d_%d_%d {", bdata->start, bdata->end, s, g);
                // graphno[s][b]: number of nodes in graph.
                for(int nd=0;nd<unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count();nd++) {
                    fprintf(stderr,"%d[start=%d end=%d];",nd,unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->start,unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->end);
                    // exon_tmp.clear();
                    // exon_tmp.push_back(no2gnode[s][g][nd]->start);
                    // exon_tmp.push_back(no2gnode[s][g][nd]->end);
                    // // exon_tmp.push_back(nd*3);
                    // // exon_tmp.push_back(nd*3+1);
                    // exonIntervals.push_back(exon_tmp);
                    // for (int i = no2gnode[s][g][nd]->start; i < no2gnode[s][g][nd]->end; i++) {
                    // 	fprintf(node_cov_bed, "chr22\t%d\t%d\tNODE\t%f\t%s\n", i, i+1, no2gnode[s][g][nd]->cov, strand_symbol.chars());
                    // }
                    int node_start = 0;
                    int node_end = 0;
                    GStr node_nd(nd);
                    GStr node_name = "Node_" + bundle_start + "_" + bundle_end + "_" + node_g + "_" + node_nd;
                    fprintf(stderr, "node_name: %s\n", node_name.chars());

                    if (nd == 0) {
                        node_start = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][1]->start-50;
                        node_end = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][1]->start;
                    } else if (nd == unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count()-1){
                        node_start = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count()-2]->end-1;
                        node_end = 	unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count()-2]->end+50;
                    } else {
                        node_start = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->start-1;
                        node_end = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->end;		
                    }


                    if (nd == 0) {
                        if(s == 0) {
                            fprintf(node_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                        } else if (s == 1) {
                            fprintf(node_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                        }
                    } else if (nd == unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count()-1){
                        if(s == 0) {
                        // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\tNODE\t%f\t+\n", no2gnode[s][g][no2gnode[s][g].Count()-2]->end, no2gnode[s][g][no2gnode[s][g].Count()-2]->end+200, 0);

                            fprintf(node_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                        } else if (s == 1) {
                            // fprintf(node_cov_pos_bed, "chr22\t%d\t%d\tNODE\t%f\t-\n", no2gnode[s][g][no2gnode[s][g].Count()-2]->end, no2gnode[s][g][no2gnode[s][g].Count()-2]->end+200, 0);

                            fprintf(node_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                        }
                    } else {
                        if(s == 0) {
                            // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\t%f\t%s\n", no2gnode[s][g][nd]->start, no2gnode[s][g][nd]->end, no2gnode[s][g][nd]->cov, strand_symbol.chars());
                            fprintf(node_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t%s\n", node_start, node_end, node_name.chars(), strand_symbol.chars());
                        } else if (s == 1) {
                            fprintf(node_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t%s\n", node_start, node_end, node_name.chars(), strand_symbol.chars());
                        }
                    }
                }

                // for(int nd=0;nd<unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count()-1;nd++) {
                // 	// fprintf(stderr,"Node %d with parents:",i);
                // 	GStr node_parent_nd(nd);
                    
                // 	for(int c=0;c<unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->child.Count();c++) {
                // 		GStr node_child_nd(unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->child[c]->nodeid);
                // 		fprintf(stderr,"%d->",nd);			
                // 		fprintf(stderr,"%d;",unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->child[c]->nodeid);
                // 		GStr junction_name = "Junc_" + bundle_start + "_" + bundle_end + "_" + node_g + "_" + node_parent_nd + "->" + node_child_nd;
                // 		fprintf(stderr, "junction_name: %s\n", junction_name.chars());
                        
                // 		int junc_start = 0;
                // 		int junc_end = 0;
                // 		if (nd == 0) {
                // 			// It's the source node.
                // 			junc_start = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][1]->start;
                // 		} else {
                // 			junc_start = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->end;
                // 		}
                // 		if (unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][ unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->child[c]->nodeid ] -> start == 0) {
                // 			// The node goes to the sink.
                // 			junc_end = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count()-2]->end;
                // 		} else {
                // 			junc_end = unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][ unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1][nd]->child[c]->nodeid ] -> start;
                // 		}
                // 		if(s == 0) {
                // 			fprintf(edge_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%d\t%s\n", junc_start, junc_end, junction_name.chars(), 10, strand_symbol.chars());
                // 		} else if (s == 1) {
                // 			fprintf(edge_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t%d\t%s\n", junc_start, junc_end, junction_name.chars(), 10, strand_symbol.chars());
                // 		}
                // 	}
                // }
                // fprintf(stderr,"}\n");
                // unispg_gp->current_gidx[s]-1 += 1;




                // sink->parent.Add(node->nodeid); // add node to sink's parents

                current_gidx[s] += 1;
            }
        }



    } else {
        for(int lclg_idx=0; lclg_idx<lclg_limit; lclg_idx++) {

            fprintf(stderr, "************************************\n");
            fprintf(stderr, "*********** Iterating lclg_idx: %d (lclg_limit: %d)********\n", lclg_idx, lclg_limit);
            fprintf(stderr, "************************************\n");


            if(no2gnode[lclg_idx].Count() == 0) {
                fprintf(stderr, "First. graph node num is 0. Pass.\n");
                continue;
            }
            // bool new_unispg = true;
            // if (new_unispg) {
            //     // Create new one.
            // } else {
            //     // Add to the original one.
            // }

            // bool more_unispg = true;
            bool more_lclg = true;
            int process_unispg = false;

            int lclg_idx_start = lclg_idx;
            int lclg_idx_end = lclg_idx;
            int unispg_idx_start = current_gidx[s];
            int unispg_idx_end = current_gidx[s];
            fprintf(stderr, "Strand: %d, More than 1 graph: %d\n", s, current_gidx[s]);

            // Here, I need to find out how many unispg & lclg should be merged into 1 new unispg.
            // while(more_unispg || more_lclg) {



	// OUT_OF_RANGE=0,
	// LASTG_COUNT_0,
	// LASTG_COUNT_N_0,
	// N_LASTG_COUNT_N_0

            LCLG_ITR_STATUS lclg_itr_status = N_LASTG_COUNT_N_0;
            while(more_lclg) {
                fprintf(stderr, "more_lclg: %d.\n", more_lclg);
                //  End the loop if lclg reaches the end.
                if (lclg_idx > lclg_limit-1) {
                    fprintf(stderr, "lclg_idx (%d, limit: %d) bigger than the limit.\n", lclg_idx, lclg_limit);
                    lclg_itr_status = OUT_OF_RANGE;
                    break;
                } else if (lclg_idx == lclg_limit-1) {
                    fprintf(stderr, ">> lclg_idx (%d, limit: %d) reach the limit.\n", lclg_idx, lclg_limit);
                    more_lclg = false;
                    process_unispg = true;
                    if((no2gnode[lclg_idx].Count() == 0)) {
                        fprintf(stderr, ">> graph node num is 0 && this is the last lclg. Process the unispg\n");
                        // The lclg graph count is 0 && this is the last lclg.
                        lclg_itr_status = LASTG_COUNT_0;
                    } else {
                        // The lclg graph count is not 0 && this is the last lclg.
                        lclg_itr_status = LASTG_COUNT_N_0;
                        // Include the last node to process.
                        // lclg_idx_end += 1;
                    }

                } else if (lclg_idx < lclg_limit-1) {
                    fprintf(stderr, ">> lclg_idx (%d, limit: %d) not yet reach the limit.\n", lclg_idx, lclg_limit);
                    more_lclg = true;
                    process_unispg = false;
                    if((no2gnode[lclg_idx].Count() == 0)) {
                        // Skip the graph
                        fprintf(stderr, ">> graph node num is 0. Go to the next lclg.\n");
                        lclg_idx_end += 1;
                        lclg_idx += 1;

                        // lclg_idx += 1;
                        // lclg_idx_start = lclg_idx;
                        // lclg_idx_end = lclg_idx;
                        continue;
                    }
                    // lclg_idx_end = lclg_idx;
                    // The lclg graph count is not 0 && this is not the last lclg.
                    lclg_itr_status = N_LASTG_COUNT_N_0;
                }

                bool sep_process_last_lclg = false;
                if (lclg_itr_status == N_LASTG_COUNT_N_0 || lclg_itr_status == LASTG_COUNT_N_0) {
                    /**********************
                    ** Printing boundaries
                    ***********************/

                    fprintf(stderr, "lclg_idx_start: %d, lclg_idx_end: %d, unispg_idx_start: %d, unispg_idx_end: %d, lclg_limit: %d\n", lclg_idx_start, lclg_idx_end, unispg_idx_start, unispg_idx_end, lclg_limit);
                    fprintf(stderr, "no2gnode: %p\n", no2gnode);
                    fprintf(stderr, "count: %d \n", no2gnode[lclg_idx].Count());
                    uint lclg_start = no2gnode[lclg_idx][1]->start;
                    fprintf(stderr, "boundary start: %d \n", lclg_start);
                    uint lclg_end = no2gnode[lclg_idx][ no2gnode[lclg_idx].Count()-2 ]->end;
                    // no2gnode->Get(no2gnode->Count()-2)->end;
                    fprintf(stderr, "boundary end: %d \n", lclg_end);
                    // while (no2gnode_unispg[s][ current_gidx[s] ].Count() == 0) {
                    //     current_gidx[s] += 1;
                    // }
                    fprintf(stderr, "$$$      unispg_idx: %d\n", current_gidx[s]);
                    fprintf(stderr, "$$$      no2gnode_unispg[s]: %d\n", sizeof(no2gnode_unispg[s]));
                    fprintf(stderr, "$$$      current_gidx[s]].Count(): %d\n", no2gnode_unispg[s][ current_gidx[s] ].Count());
                    fprintf(stderr, "$$$      no2gnode_unispg[s][current_gidx[s]][1]->start: %d\n", no2gnode_unispg[s][current_gidx[s]][1]->start);
                    uint unispg_start = no2gnode_unispg[s][current_gidx[s]][1]->start;
                    uint unispg_end = no2gnode_unispg[s][current_gidx[s]][ no2gnode_unispg[s][current_gidx[s]].Count()-2 ]->end;
                    fprintf(stderr, "$$$ lclg_start: %u,  lclg_end: %u,  unispg_start: %u,  unispg_end: %u\n", lclg_start, lclg_end, no2gnode_unispg[s][current_gidx[s]][1]->start, no2gnode_unispg[s][current_gidx[s]][ no2gnode_unispg[s][current_gidx[s]].Count()-2 ]->end);
                    for (int i = 0; i < no2gnode[lclg_idx].Count(); i++) {
                        // no2gnode->Get(i);
                        fprintf(stderr, "Before ~~ &&&& This is the local graphnode: %d\n", no2gnode[lclg_idx][i]->nodeid);
                    }


                    // unispg: -------
                    // lclg: ........
                    if (unispg_end < lclg_start) {
                        // ----------   |(s).................(e)|
                        fprintf(stderr,"\n  &&& Graph: ----------   |(s).................(e)|\n");
                        // This is the end of the lclg & unispg comparison. Only need to use the unispg. 
                        // more_unispg = true;
                        // more_lclg = false;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_unispg = true; // Process when it is the end
                        if (lclg_itr_status == LASTG_COUNT_N_0) {
                            sep_process_last_lclg = true;
                        }
                    } else if (unispg_end == lclg_start) {
                        // ----------|(s).................(e)|
                        fprintf(stderr,"\n  &&& Graph: ----------|(s).................(e)| \n");
                        // more_unispg = true;
                        // more_lclg = false;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            lclg_idx_end += 1;
                            lclg_idx += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end < lclg_end) {
                        // -----|(s)-----............(e)|
                        fprintf(stderr,"\n  &&& Graph: -----|(s)-----............(e)| \n");
                        // more_unispg = true;
                        // more_lclg = false;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            lclg_idx_end += 1;
                            lclg_idx += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end == lclg_end) {
                        // -----|(s)--------------(e)|
                        fprintf(stderr,"\n  &&& Graph: -----|(s)--------------(e)| \n");
                        //  This is the end of the lclg & unispg comparison. Set false first, but I need to check whether the next unispg and lclg overlap with the current graph.
                        // more_unispg = true;
                        // more_lclg = true;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_unispg = true; // Process when it is the end
                    } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end > lclg_end) {
                        // -----|(s)------------(e)|--
                        fprintf(stderr,"\n &&& Graph: -----|(s)------------(e)|-- \n");
                        // more_unispg = false;
                        // more_lclg = true;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
                        // |(s)----------.................(e)| 
                        fprintf(stderr,"\n &&& Graph: |(s)----------............(e)|\n");
                        // more_unispg = true;
                        // more_lclg = false;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            lclg_idx_end += 1;
                            lclg_idx += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
                        // |(s)----------(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)----------(e)| \n");
                        // This is the end of the lclg & unispg comparison. 
                        // more_unispg = true;
                        // more_lclg = true;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_unispg = true; // Process when it is the end
                    } else if (unispg_start == lclg_start && unispg_end > lclg_end) {
                        // |(s)----------(e)|-----
                        fprintf(stderr,"\n &&& Graph: |(s)----------(e)|----- \n");
                        // more_unispg = false;
                        // more_lclg = true;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end < lclg_end) {
                        // |(s)........----------........(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)........----------........(e)| \n");
                        // more_unispg = true;
                        // more_lclg = false;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            lclg_idx_end += 1;
                            lclg_idx += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end == lclg_end) {
                        // |(s)............----------(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)............----------(e)| \n");
                        // This is the end of the lclg & unispg comparison. 
                        // more_unispg = true;
                        // more_lclg = true;
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_unispg = true; // Process when it is the end
                    } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end > lclg_end) {
                        // |(s)...............------(e)|-----
                        fprintf(stderr,"\n &&& Graph: |(s)...............------(e)|-----  \n");
                        // more_unispg = false;
                        // more_lclg = true;
                        lclg_idx_end += 1; 
                        lclg_idx += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start == lclg_end) {
                        // The node is outside the current bundle => This node belongs to the next bundlenode
                        // |(s).................(e)|----------
                        fprintf(stderr,"\n &&& Graph: |(s).................(e)|----------  \n");
                        // more_unispg = false;
                        // more_lclg = true;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_unispg = false;
                        if (!more_lclg) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_unispg = true;
                        }
                    } else if (unispg_start > lclg_end) {
                        // The node is outside the current bundle => This node belongs to the next bundlenode
                        // |(s).................(e)|   ----------
                        fprintf(stderr,"\n &&& Graph: |(s).................(e)|   ---------- \n");
                        // This is the end of the lclg & unispg comparison. Only need to use the lclg.
                        // more_unispg = false;
                        // more_lclg = true;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        // [unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                        process_unispg = true; // Process when it is the end
                    } else {
                        fprintf(stderr,"\n &&& Unknown area!!!! \n");
                        process_unispg = false;
                    }

                    // If there are no lclg => process the new unispg
                    fprintf(stderr, "process_unispg: %d, more_lclg: %d\n", process_unispg, more_lclg);
                    if (!process_unispg) {
                        if (!more_lclg) {
                            fprintf(stderr, "process_unispg: %d, more_lclg: %d\n", process_unispg, more_lclg);
                            process_unispg = true;
                        }
                    }
                } else if (lclg_itr_status == LASTG_COUNT_0) {
                    lclg_idx_end += 1;
                }


                /*****************************
                 * Create a new graph 
                 *****************************/
                if (process_unispg) {

                    fprintf(stderr, "************************************\n");
                    fprintf(stderr, "*********** New UNISPG!!!!! ********\n");
                    fprintf(stderr, "************************************\n");


                    fprintf(stderr, "Inside `process_unispg`!!!\n");
                    for (int i = unispg_idx_start; i < unispg_idx_end; i++) {
                        fprintf(stderr, "unispg_index -> i: %d\n", i);
                        for (int j = 0; j < no2gnode_unispg[s][i].Count(); j++) {
                            fprintf(stderr, "** no2gnode_unispg[i][j]->start: %d\n", no2gnode_unispg[s][i][j]->start);
                            fprintf(stderr, "** no2gnode_unispg[i][j]->end: %d\n", no2gnode_unispg[s][i][j]->end);
                        }
                    }
                    for (int i = lclg_idx_start; i < lclg_idx_end; i++) {
                        fprintf(stderr, "lclg_index -> i: %d\n", i);
                        for (int j = 0; j < no2gnode[i].Count(); j++) {
                            fprintf(stderr, "** no2gnode[i][j]->start: %d\n", no2gnode[i][j]->start);
                            fprintf(stderr, "** no2gnode[i][j]->end: %d\n", no2gnode[i][j]->end);
                        }
                    }
                    if (sep_process_last_lclg) {
                        fprintf(stderr, "************************************\n");
                        fprintf(stderr, "*********** New UNISPG!!!!! ********\n");
                        fprintf(stderr, "************************************\n");
                        fprintf(stderr, "Seperately process the last lclg \n");
                        fprintf(stderr, "lclg_index -> i: %d\n", lclg_idx_end);
                        for (int j = 0; j < no2gnode[lclg_idx_end].Count(); j++) {
                            fprintf(stderr, "** no2gnode[i][j]->start: %d\n", no2gnode[lclg_idx_end][j]->start);
                            fprintf(stderr, "** no2gnode[i][j]->end: %d\n", no2gnode[lclg_idx_end][j]->end);
                        }
                    }
                    
                    fprintf(stderr, "\n\n\n");

                    unispg_idx_start = current_gidx[s];
                    unispg_idx_end = current_gidx[s];

                    // lclg_idx = lclg_idx_end;
                    lclg_idx_start = lclg_idx;
                    lclg_idx_end = lclg_idx;











                    // Creating boundary list!
                    GVec<uint> boundaries;
                    GVec<CGraphBoundaryType> boundaries_types;

                    unsigned int new_unispg_nodeid = 1;
                    // no2gnode_unispg[s][unispg_idx];
                    GPVec<CGraphnodeUnispg>* new_no2gnode_unispg; // for each graph g, on a strand s, no2gnode[g][i] gives the node i
                    new_no2gnode_unispg = new GPVec<CGraphnodeUnispg>[20000];

                    int lclg_idx = 1;
                    int unispg_idx = 1;

                    bool more_comparison = true;
                    bool lclg_is_end = false;
                    bool unispg_is_end = false;

                    GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                    GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                    GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                    CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, is_passed_s, cov_s, capacity_s, false, 0, 0, 0);

                    new_no2gnode_unispg[unispg_idx].Add(source);

                    // Iterating boundaries.
                    // int cur_boundary_start = 0;
                    // int cur_boundary_end = 0;
                    // int prev_boundary_start = 0;
                    int prev_boundary = 0;

                    CGraphnode* lclg_node =  no2gnode->Get(1);
                    CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][unispg_idx].Get(1);
                    prev_boundary = (unispg_node->start < lclg_node->start) ? unispg_node->start:lclg_node->start;
                    while(!lclg_is_end || !unispg_is_end) {
                        lclg_is_end = (lclg_idx == no2gnode->Count()-2);
                        unispg_is_end = (unispg_idx == no2gnode_unispg[s][unispg_idx].Count()-2);

                        CGraphnode* lclg_node =  no2gnode->Get(lclg_idx);
                        // fprintf(stderr, "&& lclg_node: %d (%d - %d)\t", lclg_node->nodeid, lclg_node->start, lclg_node->end);
                        CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][unispg_idx].Get(unispg_idx);
                        // fprintf(stderr, "&& unispg_node: %d (%d - %d)\t", unispg_node->nodeid, unispg_node->start, unispg_node->end);
                        CGraphnodeUnispg* node;

                        fprintf(stderr, "****** >> prev_boundary: %d \n", prev_boundary);
                        fprintf(stderr, "****** >> lclg_idx: %d,  unispg_idx: %d\n", lclg_idx, unispg_idx);
                        fprintf(stderr, "****** >> lclg_node->start: %d,  lclg_node->end: %d\n", lclg_node->start, lclg_node->end);
                        fprintf(stderr, "****** >> unispg_node->start: %d,  unispg_node->end: %d\n", unispg_node->start, unispg_node->end);

                        bool lclg_move = false;
                        bool unispg_move = false;

                        if (unispg_node->start < lclg_node->start) {
                            AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                            uint node_start_pos = 0;
                            if (prev_boundary < unispg_node->start) {
                                node_start_pos = unispg_node->start;
                            } else if (prev_boundary == unispg_node->start) {
                                node_start_pos = unispg_node->start;
                            } else if (prev_boundary > unispg_node->start) {
                                node_start_pos = prev_boundary;
                            }
                            if (unispg_node->end < lclg_node->start) {
                                fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
                                AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                prev_boundary = unispg_node->end;
                                MoveUnispg(unispg_is_end, unispg_move);
                                // if (unispg_is_end) {
                                //     AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                                //     AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                // }
                            } else if (unispg_node->end == lclg_node->start) {
                                fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
                                AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_S);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                prev_boundary = unispg_node->end;
                                MoveUnispg(unispg_is_end, unispg_move);
                            } else if (lclg_node->start < unispg_node->end) {
                                AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->start, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                if (unispg_node->end < lclg_node->end) {
                                    fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
                                    AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;

                                    prev_boundary = unispg_node->end;
                                    MoveUnispg(unispg_is_end, unispg_move);
                                } else if (unispg_node->end == lclg_node->end) {
                                    fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----------(e)|\n");
                                    AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;


                                    prev_boundary = unispg_node->end;
                                    MoveUnispg(unispg_is_end, unispg_move);
                                    MoveLclg(lclg_is_end, lclg_move);
                                } else if (lclg_node->end < unispg_node->end) {
                                    fprintf(stderr,"\n  >>>>  Graph node: -----|(s)----------(e)|-----\n");
                                    AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;

                                    prev_boundary = lclg_node->end;
                                    MoveLclg(lclg_is_end, lclg_move);
                                }
                            }
                        } else if (unispg_node->start == lclg_node->start) {
                            AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S_LCLG_S);
                            uint node_start_pos = 0;
                            if (prev_boundary < unispg_node->start) {
                                node_start_pos = unispg_node->start;
                            } else if (prev_boundary == unispg_node->start) {
                                node_start_pos = unispg_node->start;
                            } else if (prev_boundary > unispg_node->start) {
                                node_start_pos = prev_boundary;
                            }


                            if (unispg_node->end < lclg_node->end) {
                                fprintf(stderr,"\n  >>>>  Graph node: |(s)----------.......(e)|\n");
                                AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                prev_boundary = unispg_node->end;
                                MoveUnispg(unispg_is_end, unispg_move);
                            } else if (unispg_node->end == lclg_node->end) {
                                fprintf(stderr,"\n  >>>>  Graph node: |(s)------------(e)|\n");
                                AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                prev_boundary = unispg_node->end;
                                MoveUnispg(unispg_is_end, unispg_move);
                                MoveLclg(lclg_is_end, lclg_move);
                            } else if (lclg_node->end < unispg_node->end) {
                                fprintf(stderr,"\n  >>>>  Graph node: |(s)---------(e)|---\n");
                                AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                prev_boundary = lclg_node->end;
                                MoveLclg(lclg_is_end, lclg_move);
                            }
                        } else if (lclg_node->start < unispg_node->start) {
                            AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                            uint node_start_pos = 0;
                            if (prev_boundary < lclg_node->start) {
                                node_start_pos = lclg_node->start;
                            } else if (prev_boundary == lclg_node->start) {
                                node_start_pos = lclg_node->start;
                            } else if (prev_boundary > lclg_node->start) {
                                node_start_pos = prev_boundary;
                            }
                            if (unispg_node->start < lclg_node->end) {
                                AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                // Add new node
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->start, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;

                                if (unispg_node->end < lclg_node->end) {
                                    fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------....(e)|\n");
                                    AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;

                                    prev_boundary = unispg_node->end;
                                    MoveUnispg(unispg_is_end, unispg_move);
                                } else if (unispg_node->end == lclg_node->end) {
                                    fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|\n");
                                    AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;

                                    prev_boundary = unispg_node->end;
                                    MoveUnispg(unispg_is_end, unispg_move);
                                    MoveLclg(lclg_is_end, lclg_move);
                                } else if (lclg_node->end < unispg_node->end) {
                                    fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
                                    AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;
                                    
                                    prev_boundary = lclg_node->end;
                                    MoveLclg(lclg_is_end, lclg_move);
                                }
                            } else if (lclg_node->end == unispg_node->start) {
                                fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|----------\n");
                                AddBoundary(boundaries, lclg_node->end, boundaries_types, UNISPG_S_LCLG_E);
                                uint node_start_pos = 0;
                                bool create_node = true;
                                if (prev_boundary < lclg_node->start) {
                                    node_start_pos = lclg_node->start;
                                } else if (prev_boundary == lclg_node->start) {
                                    node_start_pos = lclg_node->start;
                                } else if (prev_boundary > lclg_node->start) {
                                    if (prev_boundary < lclg_node->end) {
                                        node_start_pos = prev_boundary;
                                    } else if (prev_boundary == lclg_node->end) {
                                        // Do not need to create node
                                        create_node = false;
                                    } else if (prev_boundary > lclg_node->end) {
                                        // This is an impossible case. Insane check
                                    }
                                }
                                if (create_node) {
                                    GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                    GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                    GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;
                                }

                                prev_boundary = lclg_node->end;
                                MoveLclg(lclg_is_end, lclg_move);
                            } else if (lclg_node->end < unispg_node->start) {
                                fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n");
                                AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                uint node_start_pos = 0;
                                bool create_node = true;
                                if (prev_boundary < lclg_node->start) {
                                    node_start_pos = lclg_node->start;
                                } else if (prev_boundary == lclg_node->start) {
                                    node_start_pos = lclg_node->start;
                                } else if (prev_boundary > lclg_node->start) {
                                    if (prev_boundary < lclg_node->end) {
                                        node_start_pos = prev_boundary;
                                    } else if (prev_boundary == lclg_node->end) {
                                        // Do not need to create node
                                        create_node = false;
                                    } else if (prev_boundary > lclg_node->end) {
                                        // This is an impossible case. Insane check
                                    }
                                }
                                if (create_node) {
                                    GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                    GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                    GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                    // Add new node
                                    node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                    new_no2gnode_unispg[unispg_idx].Add(node);
                                    new_unispg_nodeid += 1;
                                }

                                prev_boundary = lclg_node->end;
                                MoveLclg(lclg_is_end, lclg_move);
                                // if (lclg_is_end) {
                                //     AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                //     AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                // }
                            }
                        }

                        if (lclg_move) {
                            fprintf(stderr, "lclg_idx Move ! %d\n", lclg_idx);
                            lclg_idx += 1;
                        }
                        if (unispg_move) {
                            fprintf(stderr, "unispg_idx Move ! %d\n", unispg_idx);
                            unispg_idx += 1;
                        }
                        if (!lclg_move && !unispg_move) {
                            // Both lclg and unispg can not move forward.
                            fprintf(stderr, ">>>> Local & global graph node cannot move\n");
                            // fprintf(stderr, ">>>> Local & global graph node cannot move\n")
                            // while(!lclg_is_end || !unispg_is_end) {
                            //     lclg_is_end = (lclg_idx == no2gnode->Count()-2);
                            // }

                            if (lclg_is_end && !unispg_is_end) {
                                AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                uint node_start_pos = 0;
                                if (prev_boundary < lclg_node->end) {
                                    // it is impossible. Just the insane check
                                    GError("Wrong boundaries. Check the code!!");
                                } else if (prev_boundary == lclg_node->end) {
                                    if (unispg_node->start < prev_boundary) {
    fprintf(stderr,"\n  >>>>  Graph node: -----|(s)----------(e)|-----\n");
    fprintf(stderr,"\n  >>>>  Graph node: |(s)---------(e)|---\n");
    fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
                                        node_start_pos = prev_boundary;
                                    } else if (unispg_node->start >= prev_boundary) {
    fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|----------\n");
    fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n");
                                        node_start_pos = unispg_node->start;
                                    }
                                } else if (prev_boundary > lclg_node->end) {
                                    // This situation has happend at least one. Insert the new unispg.
                                    node_start_pos = unispg_node->start;
                                }
                                // Add new node
                                fprintf(stderr, "^^^^ node_start_pos: %d,  unispg_node->end: %d \n", node_start_pos, unispg_node->end);
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;
                                // !!!!!! Do not update the prev_boundary!!!
                                // prev_boundary = unispg_node->end;
                                // Move to next node of the unispg
                                unispg_idx += 1;
                            } else if (!lclg_is_end && unispg_is_end) {
                                AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                                AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);

                                uint node_start_pos = 0;
                                if (prev_boundary < unispg_node->end) {
                                    // it is impossible. Just the insane check
                                    GError("Wrong boundaries. Check the code!!");
                                } else if (prev_boundary == unispg_node->end) {
                                    if (lclg_node->start < prev_boundary) {
    fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
    fprintf(stderr,"\n  >>>>  Graph node: |(s)----------.......(e)|\n");
    fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------....(e)|\n");
                                        node_start_pos = prev_boundary;
                                    } else if (lclg_node->start >= prev_boundary) {
    fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
    fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
                                        node_start_pos = lclg_node->start;
                                    }
                                } else if (prev_boundary > unispg_node->end) {
                                    // This situation has happend at least one. Insert the new unispg.
                                    node_start_pos = lclg_node->start;
                                }
                                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                // Add new node
                                fprintf(stderr, "^^^^ node_start_pos: %d,  unispg_node->end: %d \n", node_start_pos, lclg_node->end);
                                node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                new_no2gnode_unispg[unispg_idx].Add(node);
                                new_unispg_nodeid += 1;
                                // !!!!!! Do not update the prev_boundary!!!
                                // prev_boundary = lclg_node->end;
                                // Move to the next node of the lclg
                                lclg_idx += 1;
                            } else if (!lclg_is_end && !unispg_is_end) {
                                // This is an insane check. If any of unispg or lclg is not end, at least 1 must move. 
                                    GError("Wrong boundaries. Check the code!!");
                            } else if (lclg_is_end && unispg_is_end) {
                                if (unispg_node->start < lclg_node->start) {
                                    if (unispg_node->end < lclg_node->start) {
                                        fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
                                        AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                                        AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                        
                                        GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                        GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                        GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                        // Add new node
                                        node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                        new_no2gnode_unispg[unispg_idx].Add(node);
                                        new_unispg_nodeid += 1;
                                        prev_boundary = lclg_node->end;

                                    } else if (unispg_node->end == lclg_node->start) {
                                        fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
                                        AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);

                                        GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                        GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                        GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                        // Add new node
                                        node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                        new_no2gnode_unispg[unispg_idx].Add(node);
                                        new_unispg_nodeid += 1;
                                        prev_boundary = lclg_node->end;

                                    } else if (lclg_node->start < unispg_node->end) {
                                        if (unispg_node->end < lclg_node->end) {
                                            fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
                                            AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                            GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                            GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                            GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                            new_no2gnode_unispg[unispg_idx].Add(node);
                                            new_unispg_nodeid += 1;
                                            prev_boundary = lclg_node->end;

                                        } else if (unispg_node->end == lclg_node->end) {
                                            fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----------(e)|\n");
                                            // Do not need to create a new node
                                        } else if (lclg_node->end < unispg_node->end) {
                                            fprintf(stderr,"\n  >>>>  Graph node: -----|(s)----------(e)|-----\n");
                                            AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                            new_no2gnode_unispg[unispg_idx].Add(node);
                                            new_unispg_nodeid += 1;
                                            prev_boundary = unispg_node->end;   
                                        }
                                    }
                                } else if (unispg_node->start == lclg_node->start) {
                                    if (unispg_node->end < lclg_node->end) {
                                        fprintf(stderr,"\n  >>>>  Graph node: |(s)----------.......(e)|\n");
                                        AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                        GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                        GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                        GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                        // Add new node
                                        node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                        new_no2gnode_unispg[unispg_idx].Add(node);
                                        new_unispg_nodeid += 1;
                                        prev_boundary = lclg_node->end;

                                    } else if (unispg_node->end == lclg_node->end) {
                                        fprintf(stderr,"\n  >>>>  Graph node: |(s)------------(e)|\n");
                                    } else if (lclg_node->end < unispg_node->end) {
                                        fprintf(stderr,"\n  >>>>  Graph node: |(s)---------(e)|---\n");
                                        AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                        // Add new node
                                        node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                        new_no2gnode_unispg[unispg_idx].Add(node);
                                        new_unispg_nodeid += 1;
                                        prev_boundary = unispg_node->end;   
                                    }
                                } else if (lclg_node->start < unispg_node->start) {
                                    if (unispg_node->start < lclg_node->end) {
                                        if (unispg_node->end < lclg_node->end) {
                                            fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------....(e)|\n");
                                            AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                            GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                            GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                            GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                            new_no2gnode_unispg[unispg_idx].Add(node);
                                            new_unispg_nodeid += 1;
                                            prev_boundary = lclg_node->end;

                                        } else if (unispg_node->end == lclg_node->end) {
                                            fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|\n");;
                                        } else if (lclg_node->end < unispg_node->end) {
                                            fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
                                            AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                            new_no2gnode_unispg[unispg_idx].Add(node);
                                            new_unispg_nodeid += 1;
                                            prev_boundary = unispg_node->end;   
                                        }
                                    } else if (lclg_node->end == unispg_node->start) {
                                        fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|----------\n");
                                        AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                        // Add new node
                                        node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                        new_no2gnode_unispg[unispg_idx].Add(node);
                                        new_unispg_nodeid += 1;
                                        prev_boundary = unispg_node->end;  
                                    } else if (lclg_node->end < unispg_node->start) {
                                        fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n \n");
                                        AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                        AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                        // Add new node
                                        node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                        new_no2gnode_unispg[unispg_idx].Add(node);
                                        new_unispg_nodeid += 1;
                                        prev_boundary = unispg_node->end;  
                                    }
                                }
                                // if (lclg_node->start < unispg_node->start) {
                                //     if (lclg_node->end < unispg_node->start) {
                                //         fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n");
                                //     } else if (lclg_node->end >= unispg_node->start) {
                                //         fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
                                //     }
                                // } else if (lclg_node->start == unispg_node->start) {
                                // } else if (lclg_node->start > unispg_node->start) {
                                //     if (lclg_node->start >unispg_node->end) {
                                //         fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
                                //     } else if (lclg_node->start == unispg_node->end) {
                                //         fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
                                //     } else if (lclg_node->start < unispg_node->end) {
                                //         fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
                                //     }
                                // }

                                break;
                            }
                        }
                    }
                    
                    GVec<bool>* is_passed_s_sink = new GVec<bool>(sample_num-1, false);
                    GVec<float>* cov_s_sink = new GVec<float>(sample_num-1, 0.0f);
                    GVec<float>* capacity_s_sink = new GVec<float>(sample_num-1, 0.0f);
                    CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, new_unispg_nodeid, is_passed_s_sink, cov_s_sink, capacity_s_sink, false, 0, 0, 0);

                    new_no2gnode_unispg[unispg_idx].Add(sink);
                    // CGraphnodeUnispgType = UNI_LCL_NODE;

                    no2gnode_unispg[s][unispg_idx].Clear();
                    no2gnode_unispg[s][unispg_idx] = new_no2gnode_unispg[unispg_idx];

                    fprintf(stderr, "Boundaries: ");
                    for (int i = 0; i < boundaries.Count(); i++) {
                        fprintf(stderr, "%d, ", boundaries.Get(i));
                    }














                }








                




                // // Check whether the next lclg is valid.
                // if (lclg_idx == lclg_limit) {
                //     more_lclg = false;
                // }
            }
        }
    }
}


void UnispgGp::AddBoundary(GVec<uint>& boundaries, uint boundary, GVec<CGraphBoundaryType>& boundaries_type, CGraphBoundaryType boundary_type) {
    if ((boundary > boundaries.Last() || boundaries.Count() == 0)) {
        boundaries.Add(boundary);
        boundaries_type.Add(boundary_type);
        fprintf(stderr,"\n  Boundary: %d\n", boundary);
        fprintf(stderr,"  >>>>  boundary_type: %s\n", enum_str[boundary_type]);
    }
}

void UnispgGp::MoveUnispg(bool& unispg_is_end, bool& unispg_move) {
    if (!unispg_is_end) {
        unispg_move = true;
    }
}

void UnispgGp::MoveLclg(bool& lclg_is_end, bool& lclg_move) {
    if (!lclg_is_end) {
        lclg_move = true;
    }
}

void UnispgGp::PrintGraphGp() {
    fprintf(stderr, "*********************************\n");
    fprintf(stderr, "*********** PrintGraphGp ********\n");
    fprintf(stderr, "*********************************\n");

    { // sDEBUG ONLY
        printTime(stderr);
        // for(int s=0;s<2;s++) {
        //     fprintf(stderr, "\n\tThere are %d stranded[%d] graphs\n", gpSize[s],int(2*s));
        //     // if () {

        //     // }
        //     for(int b=0;b<gpSize[s];b++) {
        //         fprintf(stderr, ">>>>>>> 1-2 gpSize[%d]: %d\n", s, gpSize[s]);
        //         fprintf(stderr, ">>>>>>> 1-2 graphnoGp[%d][%d]: %d\n", s, b, graphnoGp[s][b]);
        //         if(graphnoGp[s][b]) {
        //             GStr pat;
        //             fprintf(stderr,"\t\tGraph[%d][%d] with %d nodes and %d edges :",int(2*s),b,graphnoGp[s][b],edgenoGp[s][b]);
        //             for(int nd=1;nd<graphnoGp[s][b]-1;nd++) {
        //                 fprintf(stderr, "&&& nd: %d\n", nd);
        //                 fprintf(stderr, "&&& no2gnodeGp[s][b][nd]->start: %d\n", no2gnodeGp[s][b][nd]->start);
        //                 fprintf(stderr, "&&& no2gnodeGp[s][b][nd]->end: %d\n", no2gnodeGp[s][b][nd]->end);
        //                 fprintf(stderr," %d(%d-%d)",nd,no2gnodeGp[s][b][nd]->start,no2gnodeGp[s][b][nd]->end);
        //             }
        //         } 
        //         fprintf(stderr,"\n");
        //     }
        // }

    //     for(int s=0;s<2;s++) {
    // 	fprintf(stderr, "There are %d stranded[%d] graphs\n",bno[s],int(2*s));
    // 	// fprintf(stderr, "1 bundle[sno].Count(): %d\n", bno[s]);
    // 	for(int b=0;b<bno[s];b++) {
    // 		if(graphno[s][b]) {
    // 			GStr pat;
    // 			fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges with lastgpos=%d:",int(2*s),b,graphno[s][b],edgeno[s][b],lastgpos[s][b]);
    // 			for(int nd=1;nd<graphno[s][b]-1;nd++)
    // 				fprintf(stderr," %d(%d-%d)",nd,no2gnode[b][nd]->start,no2gnode[b][nd]->end);
    // 			fprintf(stderr,"\n");
    // 			print_pattern(tr2no[s][b],pat,graphno[s][b]);
    // 		}
    // 	}
    // }
    }
}

void UnispgGp::WriteGraphGp() {
    fprintf(stderr, "*********************************\n");
    fprintf(stderr, "*********** WriteGraphGp ********\n");
    fprintf(stderr, "*********************************\n");
    //  DOT file outut here 
    //  not capacity and rate 
    //  only edge weight
    for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
        int s=sno/2; // adjusted strand due to ignoring neutral strand
        fprintf(stderr, "no2gnode_unispg[s]->Count(): %d\n", no2gnode_unispg[s]->Count());
        for(int g=0;g<no2gnode_unispg[s]->Count();g++) {
            fprintf(stderr, "New writing out place!! Start writing out DOT file!!\n");
            fprintf(stderr,"after traverse:\n");
            // uinigraph_out
            fprintf(stderr,"strict digraph %d_%d_%d_%d {", no2gnode_unispg[s][g].Get(1)->start, no2gnode_unispg[s][g].Get(no2gnode_unispg[s][g].Count()-2)->end, s, g);
            // graphno[s][b]: number of nodes in graph.

            for(int nd=1;nd<no2gnode_unispg[s][g].Count()-1;nd++)
                fprintf(stderr,"%d[start=%d end=%d cov=%f];",nd,no2gnode_unispg[s][g][nd]->start,no2gnode_unispg[s][g][nd]->end,no2gnode_unispg[s][g][nd]->cov_s->Get(0));

            for(int nd=0;nd<no2gnode_unispg[s][g].Count();nd++) {
                // fprintf(stderr,"Node %d with parents:",i);
                for(int c=0;c<no2gnode_unispg[s][g][nd]->child.Count();c++) {
                    fprintf(stderr,"%d->",nd);			
                    fprintf(stderr,"%d;",no2gnode_unispg[s][g][nd]->child[c]->nodeid);
                }
            }

            fprintf(stderr,"}\n");
        }
    }
}


GPVec<CGraphnodeUnispg>** UnispgGp::get_no2gnodeGp () {
    return no2gnode_unispg;
}