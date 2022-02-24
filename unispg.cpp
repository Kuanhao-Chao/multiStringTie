#include "unispg.h"
#include "GBitVec.h"
#include <float.h>
#include <limits.h>

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

void UnispgGp::WriteLCLG(int fidx, int s, GPVec<CGraphnode>* no2gnode, int g) {
    /****************
     **  Writing out the visualization graph for the local graph.
    ****************/
    if(no2gnode[g].Count() == 0) {
        fprintf(stderr, "First. graph node num is 0. Pass.\n");
        return;
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
            node_start = no2gnode[g][1]->start-20;
            node_end = no2gnode[g][1]->start;
        } else if (nd == no2gnode[g].Count()-1){
            node_start = no2gnode[g][no2gnode[g].Count()-2]->end-1;
            node_end = 	no2gnode[g][no2gnode[g].Count()-2]->end+20;
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
}

// This function is only for checking the result of applying first UNISPG algorithm.
void UnispgGp::WriteNonOVP(int fidx, int s, int unispg_start_idx, int unispg_end_idx, GPVec<CGraphnodeUnispg>* lclg_nonoverlap) {
    /****************
     **  Writing out the visualization graph for the global graph.
     ****************/
    for (int g=unispg_start_idx; g<unispg_end_idx; g++) {
        for (int check_node=0; check_node<lclg_nonoverlap[g].Count(); check_node++) {
            fprintf(stderr, "%d, ", lclg_nonoverlap[g][check_node]->nodeid);
        }

        fprintf(stderr,"Traversing the universal splice graph!!!\n");
        // fprintf(stderr,"Unispg %d_%d_%d_%d {", bdata->start, bdata->end, s, g);
        // graphno[s][b]: number of nodes in graph.

        GStr strand_symbol;
        if (s == 0) {
            strand_symbol = "-";
        } else if (s == 1) {
            strand_symbol = "+";
        }

        for(int nd=0;nd<lclg_nonoverlap[g].Count();nd++) {
            fprintf(stderr,"%d[start=%d end=%d];",nd,lclg_nonoverlap[g][nd]->start,lclg_nonoverlap[g][nd]->end);
            int node_start = 0;
            int node_end = 0;
            GStr node_nd(nd);
            
            GStr unispg_start("");
            GStr unispg_end("");

            unispg_start = int(lclg_nonoverlap[g][1]->start);
            unispg_end = int(lclg_nonoverlap[g][ lclg_nonoverlap[g].Count()-2 ]->end);

            GStr node_name = "Node_" + unispg_start + "_" + unispg_end + "_" + node_nd;
            fprintf(stderr, "node_name: %s\n", node_name.chars());

            if (nd == 0) {
                node_start = lclg_nonoverlap[g][1]->start-20;
                node_end = lclg_nonoverlap[g][1]->start;
            } else if (nd == lclg_nonoverlap[g].Count()-1){
                node_start = lclg_nonoverlap[g][lclg_nonoverlap[g].Count()-2]->end-1;
                node_end = 	lclg_nonoverlap[g][lclg_nonoverlap[g].Count()-2]->end+20;
            } else {
                node_start = lclg_nonoverlap[g][nd]->start-1;
                node_end = lclg_nonoverlap[g][nd]->end;		
            }


            if (nd == 0) {
                if(s == 0) {
                    fprintf(node_cov_neg_novp_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                } else if (s == 1) {
                    fprintf(node_cov_pos_novp_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                }
            } else if (nd == lclg_nonoverlap[g].Count()-1){
                if(s == 0) {
                // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\tNODE\t%f\t+\n", no2gnode[s][g][no2gnode[s][g].Count()-2]->end, no2gnode[s][g][no2gnode[s][g].Count()-2]->end+200, 0);

                    fprintf(node_cov_neg_novp_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                } else if (s == 1) {
                    // fprintf(node_cov_pos_bed, "chr22\t%d\t%d\tNODE\t%f\t-\n", no2gnode[s][g][no2gnode[s][g].Count()-2]->end, no2gnode[s][g][no2gnode[s][g].Count()-2]->end+200, 0);

                    fprintf(node_cov_pos_novp_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                }
            } else {
                if(s == 0) {
                    // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\t%f\t%s\n", no2gnode[s][g][nd]->start, no2gnode[s][g][nd]->end, no2gnode[s][g][nd]->cov, strand_symbol.chars());
                    fprintf(node_cov_neg_novp_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t%s\n", node_start, node_end, node_name.chars(), strand_symbol.chars());
                } else if (s == 1) {
                    fprintf(node_cov_pos_novp_bed_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t%s\n", node_start, node_end, node_name.chars(), strand_symbol.chars());
                }
            }
        }
    }
}

void UnispgGp::WriteUNISPG(int fidx, int s, int unispg_start_idx, int unispg_end_idx) {
    /****************
     **  Writing out the visualization graph for the global graph.
     ****************/
    for (int g=unispg_start_idx; g<unispg_end_idx; g++) {
        for (int check_node=0; check_node<no2gnode_unispg[s][g].Count(); check_node++) {
            fprintf(stderr, "%d, ", no2gnode_unispg[s][g][check_node]->nodeid);
        }
        // fprintf(stderr, "\n");
        // for (int check_node=0; check_node<no2gnode_unispg[s][g].Count(); check_node++) {
        //     no2gnode_unispg[s][g][check_node]->nodeid = check_node;
        // }
        // for (int check_node=0; check_node<no2gnode_unispg[s][g].Count(); check_node++) {
        //     fprintf(stderr, "%d, ", no2gnode_unispg[s][g][check_node]->nodeid);
        // }
        // fprintf(stderr, "\n\n");

        // fprintf(stderr, "&& unispg_gp->current_gidx: %d\n", unispg_gp->current_gidx[s]-1);
        // GPVec<CGraphnode>** unispg_gp->no2gnode_unispg = unispg_gp->get_no2gnodeGp();
        // fprintf(stderr, "unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count(): %d \n", unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count());
        // int refstart_unispg = unispg_gp->no2gnode_unispg[s][g][1]->start;
        // int refend_unispg = unispg_gp->no2gnode_unispg[s][g][unispg_gp->no2gnode_unispg[s][g].Count()-2]->end;
        // GVec<int>* graphno_unispg = unispg_gp->get_graphnoGp();
        // GVec<int>* edgeno_unispg = unispg_gp->get_edgenoGp();

        fprintf(stderr,"Traversing the universal splice graph!!!\n");
        // fprintf(stderr,"Unispg %d_%d_%d_%d {", bdata->start, bdata->end, s, g);
        // graphno[s][b]: number of nodes in graph.

        GStr strand_symbol;
        if (s == 0) {
            strand_symbol = "-";
        } else if (s == 1) {
            strand_symbol = "+";
        }

        for(int nd=0;nd<no2gnode_unispg[s][g].Count();nd++) {
            fprintf(stderr,"%d[start=%d end=%d];",nd,no2gnode_unispg[s][g][nd]->start,no2gnode_unispg[s][g][nd]->end);
            int node_start = 0;
            int node_end = 0;
            GStr node_nd(nd);
            
            GStr unispg_start("");
            GStr unispg_end("");

            unispg_start = int(no2gnode_unispg[s][g][1]->start);
            unispg_end = int(no2gnode_unispg[s][g][ no2gnode_unispg[s][g].Count()-2 ]->end);

            GStr node_name = "Node_" + unispg_start + "_" + unispg_end + "_" + node_nd;
            fprintf(stderr, "node_name: %s\n", node_name.chars());

            if (nd == 0) {
                node_start = no2gnode_unispg[s][g][1]->start-20;
                node_end = no2gnode_unispg[s][g][1]->start;
            } else if (nd == no2gnode_unispg[s][g].Count()-1){
                node_start = no2gnode_unispg[s][g][no2gnode_unispg[s][g].Count()-2]->end-1;
                node_end = 	no2gnode_unispg[s][g][no2gnode_unispg[s][g].Count()-2]->end+20;
            } else {
                node_start = no2gnode_unispg[s][g][nd]->start-1;
                node_end = no2gnode_unispg[s][g][nd]->end;		
            }


            if (nd == 0) {
                if(s == 0) {
                    fprintf(node_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                } else if (s == 1) {
                    fprintf(node_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                }
            } else if (nd == no2gnode_unispg[s][g].Count()-1){
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
        // unispg_gp->current_gidx[s]-1 += 1;=
        // sink->parent.Add(node->nodeid); // add node to sink's parents
    }
}

void UnispgGp::MergeLCLG(int s, int sample_num, GPVec<CGraphnode>* no2gnode, int lclg_limit, uint boudleGP_start_idx, uint boudleGP_end_idx, int& new_nonolp_lclg_idx, bool write_unispg, GPVec<CGraphnodeUnispg>* lclg_nonoverlap) {
    fprintf(stderr, "\n*****************************\n");
    fprintf(stderr, "*********** AddGraph ********\n");
    fprintf(stderr, "*****************************\n");
    
    bool update_unispg_idx = false;
    int new_unispg_nodeid = 0;
    int process_graph_num = 0;

    for (int g_gp=boudleGP_start_idx; g_gp<boudleGP_end_idx; g_gp++) {
        fprintf(stderr, "g_gp: %d,  boudleGP_start_idx: %d,  boudleGP_end_idx: %d\n", g_gp, boudleGP_start_idx, boudleGP_end_idx);

        if(no2gnode[g_gp].Count() == 0) {
            fprintf(stderr, "Second. graph node num is 0. Pass.\n");
            continue;
        }
        // Only go inside this function once => update the unispg index.
        if (!update_unispg_idx) {
            update_unispg_idx = true;
            // If there are at least one graph => Insert source
            fprintf(stderr, "Adding source!!\n");
            GVec<bool>* is_passed_source = new GVec<bool>(sample_num-1, false);
            GVec<float>* cov_source = new GVec<float>(sample_num-1, 0.0f);
            GVec<float>* capacity_source = new GVec<float>(sample_num-1, 0.0f);
            // CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, 0, new_unispg_nodeid, is_passed_source, cov_source, capacity_source, false, 0, 0, 0);
            CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, new_unispg_nodeid, is_passed_source, cov_source, capacity_source, true, 0, 0, 0);
            new_unispg_nodeid += 1;
            if (write_unispg) {
                no2gnode_unispg[s][current_gidx[s]+new_nonolp_lclg_idx].Add(source);
            } else {
                fprintf(stderr, "lclg_nonoverlap[new_nonolp_lclg_idx]Add(source): \n");
                lclg_nonoverlap[new_nonolp_lclg_idx].Add(source);
            }
        }
        fprintf(stderr, "(%d) processed no2gnode[i].Count(): %d\n", g_gp, no2gnode[g_gp].Count());


// g_gp, i => old node id
// new_unispg_nodeid => new node id
// GHashMap<int, int> new2_nodehash(false); //hash of pointers
        for (int i=1; i<no2gnode[g_gp].Count()-1; i++) {
            CGraphnode* node = new CGraphnode(no2gnode[g_gp][i]);

            GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
            GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
            GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);

            // CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(sample_num, node->start, node->end, g_gp, i, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, node->cov, node->capacity, node->rate);
            CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(sample_num, node->start, node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, node->cov, node->capacity, node->rate);
            new_unispg_nodeid += 1;

    // // Linking parent and child
    // fprintf(stderr, "CGraphnode parent: ");
    // for(int p=0; p<node->parent.Count(); p++) {
    //     fprintf(stderr, "%d  ", node->parent[p]);
    //     // fprintf(stderr, "%d  ", no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->parent[p] ]->nodeid);

    //     // no2gnode_unispg[s][unispg_idx][i]->parent.Add(no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->parent[p] ]);
    // }
    // fprintf(stderr, "\nCGraphnode child: ");
    // for(int c=0; c<no2gnode[g][i]->child.Count(); c++) {
    //     fprintf(stderr, "%d  ", no2gnode[g][i]->child[c]); 
    //     // fprintf(stderr, "%d  ", no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->child[c] ]->nodeid); 

    //     // no2gnode_unispg[s][unispg_idx][i]->child.Add(no2gnode_unispg[s][unispg_idx][ no2gnode[g][i]->child[c] ]);
    // }
    // fprintf(stderr, "\n");
            if (write_unispg) {
                fprintf(stderr, "Writing out unispg\n");
                int g_idx_tmp = current_gidx[s]+new_nonolp_lclg_idx;
                if (no2gnode_unispg[s][g_idx_tmp].Count() == 1) {
                    no2gnode_unispg[s][g_idx_tmp].Add(node_unispg);
                } else {
                    fprintf(stderr, "no2gnode_unispg[s][current_gidx[s]+new_nonolp_lclg_idx].Count(): %d\n", no2gnode_unispg[s][g_idx_tmp].Count());
                    for (int insert_idx=1; insert_idx<no2gnode_unispg[s][g_idx_tmp].Count(); insert_idx++) {
                        fprintf(stderr, "(%d) node_unispg->start: %d  /  no2gnode_unispg[s][current_gidx[s]+new_nonolp_lclg_idx][insert_idx]->start: %d\n", insert_idx, node_unispg->start, no2gnode_unispg[s][g_idx_tmp][insert_idx]->start);
                        if (node_unispg->start < no2gnode_unispg[s][g_idx_tmp][insert_idx]->start) {
                            no2gnode_unispg[s][g_idx_tmp].Insert(insert_idx, node_unispg); 
                            break;
                        } else {
                            if (insert_idx == no2gnode_unispg[s][g_idx_tmp].Count()-1) {
                                no2gnode_unispg[s][g_idx_tmp].Add(node_unispg);
                                break; 
                            }
                        }
                    }
                }
            } else {
                fprintf(stderr, "Writing out unispg\n") ;
                if (lclg_nonoverlap[new_nonolp_lclg_idx].Count() == 1) {
                    fprintf(stderr, "Adding node_unispg: \n");
                    lclg_nonoverlap[new_nonolp_lclg_idx].Add(node_unispg);
                } else {
                    fprintf(stderr, "lclg_nonoverlap[new_nonolp_lclg_idx].Count(): %d\n", lclg_nonoverlap[new_nonolp_lclg_idx].Count());
                    for (int insert_idx=1; insert_idx<lclg_nonoverlap[new_nonolp_lclg_idx].Count(); insert_idx++) {
                        fprintf(stderr, "(%d) node_unispg->start: %d  /  lclg_nonoverlap[new_nonolp_lclg_idx][insert_idx].start: %d\n", insert_idx, node_unispg->start, lclg_nonoverlap[new_nonolp_lclg_idx][insert_idx]->start);
                        if (node_unispg->start < lclg_nonoverlap[new_nonolp_lclg_idx][insert_idx]->start) {
                            fprintf(stderr, "Insert node_unispg: \n");
                            lclg_nonoverlap[new_nonolp_lclg_idx].Insert(insert_idx, node_unispg); 
                            break;
                        } else {
                            if (insert_idx == lclg_nonoverlap[new_nonolp_lclg_idx].Count()-1) {
                                fprintf(stderr, "lclg_nonoverlap[new_nonolp_lclg_idx].Add(node_unispg); \n");
                                lclg_nonoverlap[new_nonolp_lclg_idx].Add(node_unispg);
                                break; 
                            }
                        }
                    }
                }
            }
        }
        // Count how many graphs are processed.
        process_graph_num += 1;
    }

    if (update_unispg_idx) {
        // Insert sink
        fprintf(stderr, "Adding sink!!\n");
        GVec<bool>* is_passed_sink = new GVec<bool>(sample_num-1, false);
        GVec<float>* cov_sink = new GVec<float>(sample_num-1, 0.0f);
        GVec<float>* capacity_sink = new GVec<float>(sample_num-1, 0.0f);
        // CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, 0, 0, new_unispg_nodeid, is_passed_sink, cov_sink, capacity_sink, false, 0, 0, 0);
        CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, new_unispg_nodeid, is_passed_sink, cov_sink, capacity_sink, true, 0, 0, 0);

        if (write_unispg) {
            no2gnode_unispg[s][current_gidx[s]+new_nonolp_lclg_idx].Add(sink);
            fprintf(stderr, "current_gidx[s]+new_nonolp_lclg_idx: %d.\n", current_gidx[s]+new_nonolp_lclg_idx);
        } else {
            lclg_nonoverlap[new_nonolp_lclg_idx].Add(sink);
        }
        new_nonolp_lclg_idx += 1;
        // current_gidx[s]+=1;
    }
}


void UnispgGp::FirstUnispgAlgo(int fidx, int s, int sample_num, GPVec<CGraphnode>* no2gnode, int lclg_limit, int& new_nonolp_lclg_idx, bool write_unispg, GPVec<CGraphnodeUnispg>* lclg_nonoverlap) {
    uint unispg_start = 0;
    uint unispg_end = 0;
    uint boudleGP_start_idx = 0;
    uint boudleGP_end_idx = 0;

    for(int g=0; g<lclg_limit; g++) {
        // fprintf(stderr, "g: %d, lclg_limit: %d\n", g, lclg_limit);
        // If the last graph is 0 => still need to process once the previous & current bundleGp
        if(no2gnode[g].Count() == 0 && g != (lclg_limit-1)) {
            // fprintf(stderr, "First. graph node num is 0. Pass.\n");
            continue;
        }
        WriteLCLG(fidx, s, no2gnode, g);
        // 1. no2gnode[g].Count() != 0
        // 2. g == (lclg_limit-1)
        boudleGP_end_idx = g;
        // Check whether the current lclg overlap with the previous lclg
        uint cur_node_start;
        uint cur_node_end;
        if(no2gnode[g].Count() == 0 && g == (lclg_limit-1)) {
            // This is the end. Process the previous group.
            cur_node_start = UINT_MAX;
            cur_node_end = UINT_MAX;
        } else {
            cur_node_start = no2gnode[g][1]->start;
            cur_node_end = no2gnode[g][ no2gnode[g].Count()-2 ]->end;
        }

        // fprintf(stderr, "g: %d, cur_node_start: %u, cur_node_end: %u, lclg_limit: %d\n", g, cur_node_start, cur_node_end, lclg_limit);
        // fprintf(stderr, "unispg_start: %d, unispg_end: %d\n", unispg_start, unispg_end);

        if (unispg_end < cur_node_start) {
            // |************|   |============|
            // Process the previous bundleGp
            MergeLCLG(s, sample_num, no2gnode, lclg_limit, boudleGP_start_idx, boudleGP_end_idx, new_nonolp_lclg_idx, write_unispg, lclg_nonoverlap);
            boudleGP_start_idx = g;
            // fprintf(stderr, "Update boudleGP_start_idx: %d.\n", boudleGP_start_idx);
            unispg_start = cur_node_start;
            unispg_end = cur_node_end;

            // Boundary case => last node
            if (g == (lclg_limit-1)) {
                if(no2gnode[g].Count() == 0) {
                    // (1) ------  xxxxxxxxxxx
                    // (2) xxxxxxxxxxxxxxxxxxx
                    // All good. No actions needed.
                } else {
                    // Need to process the last orphan node.
                    // (3) xxxxxxxxxxxxx  -------
                    // (4) ------  xxxxx  -------
                    MergeLCLG(s, sample_num, no2gnode, lclg_limit, lclg_limit-1, lclg_limit, new_nonolp_lclg_idx, write_unispg, lclg_nonoverlap);
                }
            }
            fprintf(stderr, "\nEnd\n");
        } else {
            // |************||============|
            // |********|===|========|
            // |******======|
            // |***======***|
            fprintf(stderr, "graph overlaps!!!\n");
            if (unispg_end < cur_node_end) {
                unispg_end = cur_node_end;
            }
            // unispg_end unchanged.
            // Boundary case => last node => Process the lclg_gp
            if (g == (lclg_limit-1)) {
                // Process the lclg_gp
                MergeLCLG(s, sample_num, no2gnode, lclg_limit, boudleGP_start_idx, lclg_limit, new_nonolp_lclg_idx, write_unispg, lclg_nonoverlap);
                boudleGP_start_idx = g;
                fprintf(stderr, "Update boudleGP_start_idx: %d.\n", boudleGP_start_idx);
            }
        }
    }
}




void UnispgGp::AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode, int lclg_limit) {
    int sample_num = unispg_gp->samples.Count();
    if (fidx == 0) {
        int new_nonolp_lclg_idx = 0;
        FirstUnispgAlgo(fidx, s, sample_num, no2gnode, lclg_limit, new_nonolp_lclg_idx, true, NULL);
        // WriteNonOVP(fidx, s, 0, new_nonolp_lclg_idx, no2gnode_unispg[s]);
        WriteUNISPG(fidx, s, current_gidx[s], current_gidx[s]+new_nonolp_lclg_idx);
        current_gidx[s] += new_nonolp_lclg_idx;
    } else {
        // for(int g=0; g<lclg_limit; g++) {
        //     WriteLCLG(fidx, s, no2gnode, g);
        // }
        int new_nonolp_lclg_idx = 0;
        GPVec<CGraphnodeUnispg>* lclg_nonoverlap = new GPVec<CGraphnodeUnispg>[lclg_limit];
        FirstUnispgAlgo(fidx, s, sample_num, no2gnode, lclg_limit, new_nonolp_lclg_idx, false, lclg_nonoverlap);
        WriteNonOVP(fidx, s, 0, new_nonolp_lclg_idx, lclg_nonoverlap);

        // WriteUNISPG(fidx, s, 0, new_nonolp_lclg_idx);

        // // New algorth. Assuming all graphs might be overlapping.
        // for(int lclg_idx=0; lclg_idx<lclg_limit; lclg_idx++) {
        //     fprintf(stderr, "************************************\n");
        //     fprintf(stderr, "*********** Iterating lclg_idx: %d (lclg_limit: %d)********\n", lclg_idx, lclg_limit);
        //     fprintf(stderr, "************************************\n");
        //     if(no2gnode[lclg_idx].Count() == 0) {
        //         fprintf(stderr, "First. graph node num is 0. Pass.\n");
        //         continue;
        //     }
        //     // if the current lcl graph overlap the next graph.
        // }


        
        // Old algorithm. Assuming bundles (splice graphs) do not overlap.
        for (int lclg_idx=0; lclg_idx<new_nonolp_lclg_idx; lclg_idx++) {
        // for(int lclg_idx=0; lclg_idx<lclg_limit; lclg_idx++) {

            fprintf(stderr, "************************************\n");
            fprintf(stderr, "*********** Iterating lclg_idx: %d (lclg_limit: %d)********\n", lclg_idx, new_nonolp_lclg_idx);
            fprintf(stderr, "************************************\n");
            if(lclg_nonoverlap[lclg_idx].Count() == 0) {
                fprintf(stderr, "First. graph node num is 0. Pass. (lclg_idx: %d)\n", lclg_idx);
                continue;
            }

            // bool more_unispg = true;
            bool more_lclg = true;
            bool try_more_unispg = true;
            int process_ovp_graphs = false;

            int lclg_idx_start = lclg_idx;
            int lclg_idx_end = lclg_idx;
            int unispg_idx_start = current_gidx[s];
            int unispg_idx_end = current_gidx[s];
            fprintf(stderr, "Strand: %d, More than 1 graph: %d\n", s, current_gidx[s]);

            /****************
             **  Here, I need to find out how many unispg & lclg should be merged into 1 new unispg.
            ****************/
            LCLG_ITR_STATUS lclg_itr_status = N_LASTG_COUNT_N_0;
            while(more_lclg || try_more_unispg) {
                fprintf(stderr, "\n\n^^^^^^^ Recruiting more lclg & unispg!!!!\n");
                fprintf(stderr, "^^^^^^^ more_lclg: %d,  try_more_unispg: %d\n", more_lclg, try_more_unispg);
                lclg_itr_status = LCLG_ITR_INIT;
                more_lclg = false;
                try_more_unispg = false;
                process_ovp_graphs = false;
                //  End the loop if lclg reaches the end.
                if (lclg_idx > new_nonolp_lclg_idx-1) {
                    fprintf(stderr, ">> 'OUT_OF_RANGE' lclg_idx (%d, limit: %d) bigger than the limit.\n", lclg_idx, new_nonolp_lclg_idx);
                    more_lclg = false;
                    lclg_itr_status = OUT_OF_RANGE;
                    break;
                } else if (lclg_idx == new_nonolp_lclg_idx-1) {
                    fprintf(stderr, ">> lclg_idx (%d, limit: %d) reach the limit.\n", lclg_idx, new_nonolp_lclg_idx);
                    more_lclg = false;
                    // lclg_idx_end += 1;
                    if((lclg_nonoverlap[lclg_idx].Count() == 0)) {
                        fprintf(stderr, ">> 'LASTG_COUNT_0'. graph node num is 0 && this is the last lclg. Process the unispg\n");
                        // The lclg graph count is 0 && this is the last lclg.
                        lclg_itr_status = LASTG_COUNT_0;
                    } else {
                        fprintf(stderr, ">> 'LASTG_COUNT_N_0'. graph node num is not 0 && this is the last lclg. Process the unispg\n");
                        // The lclg graph count is not 0 && this is the last lclg.
                        // Include the last node to process.
                        lclg_itr_status = LASTG_COUNT_N_0;
                    }

                } else if (lclg_idx < new_nonolp_lclg_idx-1) {
                    fprintf(stderr, ">> lclg_idx (%d, limit: %d) not yet reach the limit.\n", lclg_idx, new_nonolp_lclg_idx);
                    more_lclg = true;
                    if((lclg_nonoverlap[lclg_idx].Count() == 0)) {
                        // Skip the graph
                        fprintf(stderr, ">> 'N_LASTG_COUNT_0'. graph node num is 0. Go to the next lclg.\n");
                        lclg_itr_status = N_LASTG_COUNT_0;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        continue;
                    }
                    // The lclg graph count is not 0 && this is not the last lclg.
                    lclg_itr_status = N_LASTG_COUNT_N_0;
                }

                // bool sep_process_last_lclg = false;
                if (lclg_itr_status == N_LASTG_COUNT_N_0 || lclg_itr_status == LASTG_COUNT_N_0) {
                    // The graph is not empty.
                    uint lclg_start = lclg_nonoverlap[lclg_idx][1]->start;
                    uint lclg_end = lclg_nonoverlap[lclg_idx][ lclg_nonoverlap[lclg_idx].Count()-2 ]->end;
                    uint unispg_start = no2gnode_unispg[s][current_gidx[s]][1]->start;
                    uint unispg_end = no2gnode_unispg[s][current_gidx[s]][ no2gnode_unispg[s][current_gidx[s]].Count()-2 ]->end;

                    /**********************
                    ** Printing boundaries
                    ***********************/
                    fprintf(stderr, "** lclg_idx: %d (%d, %d), new_nonolp_lclg_idx: %d\n", lclg_idx, lclg_idx_start, lclg_idx_end, new_nonolp_lclg_idx);
                    fprintf(stderr, "** unispg_idx: %d (%d, %d), new_nonolp_lclg_idx: %d\n", current_gidx[s], unispg_idx_start, unispg_idx_end, new_nonolp_lclg_idx);
                    // fprintf(stderr, "no2gnode: %p\n", lclg_nonoverlap);
                    // fprintf(stderr, "count: %d \n", lclg_nonoverlap[lclg_idx].Count());
                    // fprintf(stderr, "boundary start: %d \n", lclg_start);
                    // // lclg_nonoverlap->Get(lclg_nonoverlap->Count()-2)->end;
                    // fprintf(stderr, "boundary end: %d \n", lclg_end);
                    // fprintf(stderr, "$$$      unispg_idx: %d\n", current_gidx[s]);
                    // fprintf(stderr, "$$$      no2gnode_unispg[s]: %d\n", sizeof(no2gnode_unispg[s]));
                    // fprintf(stderr, "$$$      current_gidx[s]].Count(): %d\n", no2gnode_unispg[s][ current_gidx[s] ].Count());
                    // fprintf(stderr, "$$$      no2gnode_unispg[s][current_gidx[s]][1]->start: %d\n", no2gnode_unispg[s][current_gidx[s]][1]->start);
                    // fprintf(stderr, "$$$ lclg_start: %u,  lclg_end: %u,  unispg_start: %u,  unispg_end: %u\n", lclg_start, lclg_end, no2gnode_unispg[s][current_gidx[s]][1]->start, no2gnode_unispg[s][current_gidx[s]][ no2gnode_unispg[s][current_gidx[s]].Count()-2 ]->end);
                    // for (int i = 0; i < lclg_nonoverlap[lclg_idx].Count(); i++) {
                    //     fprintf(stderr, "Before ~~ &&&& This is the local graphnode: %d\n", lclg_nonoverlap[lclg_idx][i]->nodeid);
                    // }


                    // unispg: -------
                    // lclg: ........
                    // When processing the graphs, 'XX_idx_end' is not included. (XX_idx < XX_idx_end)
                    if (unispg_end < lclg_start) {
                        // ----------   |(s).................(e)|
                        fprintf(stderr,"\n  &&& Graph: ----------   |(s).................(e)|\n");
                        // This is the end of the lclg & unispg comparison. Only need to use the unispg. 
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_ovp_graphs = true; // Process when it is the end
                        try_more_unispg = true;
                    } else if (unispg_end == lclg_start) {
                        // ----------|(s).................(e)|
                        fprintf(stderr,"\n  &&& Graph: ----------|(s).................(e)| \n");
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_ovp_graphs = false;
                        try_more_unispg = true;
                    } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end < lclg_end) {
                        // -----|(s)-----............(e)|
                        fprintf(stderr,"\n  &&& Graph: -----|(s)-----............(e)| \n");
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_ovp_graphs = false;
                        try_more_unispg = true;
                    } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end == lclg_end) {
                        // -----|(s)--------------(e)|
                        fprintf(stderr,"\n  &&& Graph: -----|(s)--------------(e)| \n");
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_ovp_graphs = true; // Process when it is the end
                        try_more_unispg = true;
                    } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end > lclg_end) {
                        // -----|(s)------------(e)|--
                        fprintf(stderr,"\n &&& Graph: -----|(s)------------(e)|-- \n");
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_ovp_graphs = false;
                        try_more_unispg = false;
                        if (lclg_itr_status == LASTG_COUNT_N_0) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_ovp_graphs = true;
                        }
                    } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
                        // |(s)----------.................(e)| 
                        fprintf(stderr,"\n &&& Graph: |(s)----------............(e)|\n");
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_ovp_graphs = false;
                        try_more_unispg = true;
                    } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
                        // |(s)----------(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)----------(e)| \n");
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        try_more_unispg = true;
                        process_ovp_graphs = true;
                    } else if (unispg_start == lclg_start && unispg_end > lclg_end) {
                        // |(s)----------(e)|-----
                        fprintf(stderr,"\n &&& Graph: |(s)----------(e)|----- \n");
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        process_ovp_graphs = false;
                        try_more_unispg = false;
                        if (lclg_itr_status == LASTG_COUNT_N_0) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_ovp_graphs = true;
                        }
                    } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end < lclg_end) {
                        // |(s)........----------........(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)........----------........(e)| \n");
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        process_ovp_graphs = false;
                        try_more_unispg = true;
                    } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end == lclg_end) {
                        // |(s)............----------(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)............----------(e)| \n");
                        // This is the end of the lclg & unispg comparison. 
                        unispg_idx_end += 1;
                        current_gidx[s] += 1;
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        try_more_unispg = true;
                        process_ovp_graphs = true;
                    } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end > lclg_end) {
                        // |(s)...............------(e)|-----
                        fprintf(stderr,"\n &&& Graph: |(s)...............------(e)|-----  \n");
                        lclg_idx_end += 1; 
                        lclg_idx += 1;
                        try_more_unispg = false;
                        process_ovp_graphs = false;
                        if (lclg_itr_status == LASTG_COUNT_N_0) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_ovp_graphs = true;
                        }
                    } else if (unispg_start == lclg_end) {
                        // |(s).................(e)|----------
                        fprintf(stderr,"\n &&& Graph: |(s).................(e)|----------  \n");
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        try_more_unispg = false;
                        process_ovp_graphs = false;
                        if (lclg_itr_status == LASTG_COUNT_N_0) {
                            unispg_idx_end += 1;
                            current_gidx[s] += 1;
                            process_ovp_graphs = true;
                        }
                    } else if (unispg_start > lclg_end) {
                        // The node is outside the current bundle => This node belongs to the next bundlenode
                        // |(s).................(e)|   ----------
                        fprintf(stderr,"\n &&& Graph: |(s).................(e)|   ---------- \n");
                        // This is the end of the lclg & unispg comparison. Only need to use the lclg.
                        lclg_idx_end += 1;
                        lclg_idx += 1;
                        try_more_unispg = false;
                        process_ovp_graphs = true;
                    } else {
                        fprintf(stderr,"\n &&& Unknown area!!!! \n");
                        process_ovp_graphs = false;
                    }
                } else if (lclg_itr_status == LASTG_COUNT_0) {
                    // The graph is empty.
                    process_ovp_graphs = true; // Process when it is the end
                    unispg_idx_end += 1;
                }


                /*****************************
                 * Create a new graph 
                 *****************************/
                if (process_ovp_graphs) {

                    fprintf(stderr, "************************************\n");
                    fprintf(stderr, "*********** New UNISPG!!!!! ********\n");
                    fprintf(stderr, "************************************\n");
                    /******************************************
                     * These variables are for creating CGraphnodeUnispg in the new unispg.
                     *****************************************/
                    unsigned int new_unispg_nodeid = 1;
                    GPVec<CGraphnodeUnispg>* new_no2gnode_unispg; // for each graph g, on a strand s, no2gnode[g][i] gives the node i
                    new_no2gnode_unispg = new GPVec<CGraphnodeUnispg>[20000];

                    bool more_comparison = true;

                    GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                    GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                    GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                    CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, is_passed_s, cov_s, capacity_s, false, 0, 0, 0);
                    new_no2gnode_unispg->Add(source);

                    // uint prev_boundary = 0;
                    /******************************************
                     * These variables are for iterating unispg & lclg
                     *****************************************/
                    // Graph related variable.
                    //  Tracking graph id.
                    int lclg_i=lclg_idx_start;
                    int unispg_i=unispg_idx_start;

                    uint lclg_start_pcs = 0;
                    uint lclg_end_pcs = 0;
                    uint unispg_start_pcs = 0;
                    uint unispg_end_pcs = 0;

                    // Graphnode related variable.
                    CGraphnodeUnispg* lclg_node;
                    CGraphnodeUnispg* unispg_node;
                    //  Tracking node id.
                    int lclg_node_idx = 1;
                    int unispg_node_idx = 1;
                    bool lclg_is_lastnode = false;
                    bool unispg_is_lastnode = false;


                    while (lclg_i < lclg_idx_end || unispg_i < unispg_idx_end) {
                        fprintf(stderr, "(lclg_i, lclg_idx_end) (%d, %d): %d\n", lclg_i, lclg_idx_end, lclg_i < lclg_idx_end);
                        fprintf(stderr, "(unispg_i,  unispg_idx_end) (%d, %d): %d\n",unispg_i, unispg_idx_end, unispg_i < unispg_idx_end);
                        fprintf(stderr, "&&&&&&&& lclg_i: %d,  lclg_node_idx: %d,  lclg_is_lastnode: %d\n", lclg_i, lclg_node_idx, lclg_is_lastnode);
                        fprintf(stderr, "&&&&&&&& unispg_i: %d,  unispg_node_idx: %d,  unispg_is_lastnode: %d\n", unispg_i, unispg_node_idx, unispg_is_lastnode);
                        fprintf(stderr, "&&&&&&&& prev_boundary: %d\n", prev_boundary);
                        // indicator of whether to go to the next graph. 
                        bool lclg_next = false;
                        bool unispg_next = false;

                        lclg_is_lastnode = false;
                        unispg_is_lastnode = false;
                        // fprintf(stderr, "lclg_idx_start: %d, lclg_idx_end: %d, lclg_i: %d, unispg_idx_start: %d, unispg_idx_end: %d, unispg_i: %d, lclg_limit: %d\n", lclg_idx_start, lclg_idx_end, lclg_i, unispg_idx_start, unispg_idx_end, unispg_i, lclg_limit);

                        if(lclg_nonoverlap[lclg_i].Count() == 0  && lclg_i < lclg_idx_end) {
                            fprintf(stderr, "First. graph node num is 0. Pass. (lclg_idx: %d)\n", lclg_idx);
                            lclg_next = true;
                // if (lclg_next && lclg_i < lclg_idx_end) {
                //     fprintf(stderr, "lclg_next: %d, lclg_i: %d, lclg_idx_start: %d, lclg_idx_end %d \n", lclg_next, lclg_i, lclg_idx_start, lclg_idx_end);
                //     fprintf(stderr, "unispg_next: %d, unispg_i: %d, unispg_idx_start: %d, unispg_idx_end %d \n", unispg_next, unispg_i, unispg_idx_start, unispg_idx_end);
                //     lclg_i += 1;
                //     lclg_node_idx = 1;
                // }
                        } else {
                            fprintf(stderr, "1 **** >> (%d) lclg_idx_start: %u, lclg_idx_end %u, new_nonolp_lclg_idx: %d \n", lclg_i, lclg_idx_start, lclg_idx_end, new_nonolp_lclg_idx);
                            fprintf(stderr, "1 **** >> (%d) unispg_idx_start: %u, unispg_idx_end %u \n", unispg_i, unispg_idx_start, unispg_idx_end);
                            fprintf(stderr, "1    >>> lclg boundary start: %u \n", lclg_start_pcs);
                            fprintf(stderr, "1    >>> lclg boundary end: %u \n", lclg_end_pcs);
                            fprintf(stderr, "1    >>> unispg boundary start: %u \n", unispg_start_pcs);
                            fprintf(stderr, "1    >>> unispg boundary start: %u \n", unispg_end_pcs);

                            if (lclg_i < lclg_idx_end && unispg_i < unispg_idx_end) {
                                fprintf(stderr, "   >>> Insie 'lclg_i < lclg_idx_end && unispg_i < unispg_idx_end'\n");
                                lclg_start_pcs = lclg_nonoverlap[lclg_i][1]->start;
                                lclg_end_pcs = lclg_nonoverlap[lclg_i][ lclg_nonoverlap[lclg_i].Count()-2 ]->end;
                                lclg_node = lclg_nonoverlap[lclg_i][lclg_node_idx];   
                                
                                unispg_start_pcs = no2gnode_unispg[s][unispg_i][1]->start;
                                unispg_end_pcs = no2gnode_unispg[s][unispg_i][ no2gnode_unispg[s][unispg_i].Count()-2 ]->end;
                                unispg_node = no2gnode_unispg[s][unispg_i][unispg_node_idx];


                                /******************************************
                                 * Decide to move lclg or unispg by the graph boundaries.
                                 *****************************************/
                                if (lclg_end_pcs < unispg_end_pcs) {
                                    fprintf(stderr, "&&&& Mover lclg. \n");
                                    lclg_next = true;
                                } else if (lclg_end_pcs == unispg_end_pcs) {
                                    fprintf(stderr, "&&&& Mover lclg & unispg. \n");
                                    lclg_next = true;
                                    unispg_next = true;
                                } else if (lclg_end_pcs > unispg_end_pcs) {
                                    fprintf(stderr, "&&&& Mover unispg (has_lclg && has_unispg). \n");
                                    unispg_next = true;
                                } 
                            } else if (lclg_i < lclg_idx_end && !(unispg_i < unispg_idx_end)) {
                                fprintf(stderr, "   >>> Insie 'lclg_i < lclg_idx_end && !(unispg_i < unispg_idx_end)'\n");
                                lclg_next = true;

                                lclg_start_pcs = lclg_nonoverlap[lclg_i][1]->start;
                                lclg_end_pcs = lclg_nonoverlap[lclg_i][ lclg_nonoverlap[lclg_i].Count()-2 ]->end;
                                lclg_node = lclg_nonoverlap[lclg_i][lclg_node_idx];   

                                unispg_start_pcs = UINT_MAX;
                                unispg_end_pcs = UINT_MAX;
                                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                unispg_node = new CGraphnodeUnispg(sample_num, unispg_start_pcs, unispg_end_pcs, -1, is_passed_s, cov_s, capacity_s);

                            } else if (!(lclg_i < lclg_idx_end) && unispg_i < unispg_idx_end) {
                                fprintf(stderr, "   >>> Insie '!(lclg_i < lclg_idx_end) && unispg_i < unispg_idx_end'\n");
                                // if (unispg_i+1 == unispg_idx_end) break;
                                unispg_next = true;
                                lclg_start_pcs = UINT_MAX;
                                lclg_end_pcs = UINT_MAX;
                                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                lclg_node = new CGraphnodeUnispg(sample_num, lclg_start_pcs, lclg_end_pcs, -1, is_passed_s, cov_s, capacity_s);

                                
                                unispg_start_pcs = no2gnode_unispg[s][unispg_i][1]->start;
                                unispg_end_pcs = no2gnode_unispg[s][unispg_i][ no2gnode_unispg[s][unispg_i].Count()-2 ]->end;
                                unispg_node = no2gnode_unispg[s][unispg_i][unispg_node_idx];
                            } else if (!(lclg_i < lclg_idx_end) && !(unispg_i < unispg_idx_end)) {
                                fprintf(stderr, "   >>> Insie '!(lclg_i < lclg_idx_end) && !(unispg_i < unispg_idx_end)'\n");
                                break;
                            }

                            fprintf(stderr, "2 **** >> (%d) lclg_idx_start: %u, lclg_idx_end %u, lclg_limit: %d \n", lclg_i, lclg_idx_start, lclg_idx_end, lclg_limit);
                            fprintf(stderr, "2 **** >> (%d) unispg_idx_start: %u, unispg_idx_end %u \n", unispg_i, unispg_idx_start, unispg_idx_end);
                            fprintf(stderr, "2    >>> lclg boundary start: %u \n", lclg_start_pcs);
                            fprintf(stderr, "2    >>> lclg boundary end: %u \n", lclg_end_pcs);
                            fprintf(stderr, "2    >>> unispg boundary start: %u \n", unispg_start_pcs);
                            fprintf(stderr, "2    >>> unispg boundary start: %u \n", unispg_end_pcs);

                            
                            uint node_start_pos = 0;
                            if (prev_boundary == 0) {
                                prev_boundary = (unispg_node->start < lclg_node->start) ? unispg_node->start:lclg_node->start;
                            }
                            /******************************
                            ** Iterate the nodes in lclg & unispg. If any of the graph reach the last node,
                            **  Stop the while loop and move on to the next graph pair comparison.
                            *******************************/
                            while(!lclg_is_lastnode && !unispg_is_lastnode) {
                                if (lclg_i < lclg_idx_end) {
                                    fprintf(stderr, "Node size: %d\n", lclg_nonoverlap[lclg_i].Count()-2);
                                    lclg_node = lclg_nonoverlap[lclg_i][lclg_node_idx];   
                                    lclg_is_lastnode = (lclg_node_idx == lclg_nonoverlap[lclg_i].Count()-2);
                                } 
                                else {
                                    /******************************
                                    ** Here, 'lclg_i == lclg_idx_end'
                                    **  Handle the boundary case
                                    *******************************/
                                    // lclg_is_lastnode = true;
                                }
                                if (unispg_i < unispg_idx_end) {
                                    fprintf(stderr, "Node size: %d\n", no2gnode_unispg[s][unispg_i].Count()-2);
                                    unispg_node = no2gnode_unispg[s][unispg_i][unispg_node_idx];
                                    unispg_is_lastnode = (unispg_node_idx == no2gnode_unispg[s][unispg_i].Count()-2);
                                } 
                                else {
                                    /******************************
                                    ** Here, 'unispg_i == unispg_idx_end'
                                    **  Handle the boundary case
                                    *******************************/
                                    // unispg_is_lastnode = true;
                                }
                                fprintf(stderr, "Inside iterating nodes in lclg & unispg.\n");
                                fprintf(stderr, "lclg_i: %d,  lclg_node_idx: %d,  lclg_is_lastnode: %d\n", lclg_i, lclg_node_idx, lclg_is_lastnode);
                                fprintf(stderr, "unispg_i: %d,  unispg_node_idx: %d,  unispg_is_lastnode: %d\n", unispg_i, unispg_node_idx, unispg_is_lastnode);


                                // CGraphnodeUnispg* lclg_node = lclg_nonoverlap[lclg_i][lclg_node_idx];
                                // CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][unispg_i][unispg_node_idx];
                                fprintf(stderr, "\t&& (%d) lclg_node(%d): %d,  (%u - %u)\n", lclg_i, lclg_node->nodeid, lclg_node_idx, lclg_node->start, lclg_node->end);
                                fprintf(stderr, "\t&& (%d) unispg_node(%d): %d,  (%u - %u)\n", unispg_i, unispg_node->nodeid, unispg_node_idx, unispg_node->start, unispg_node->end);

                                fprintf(stderr, "\t****** >> prev_boundary: %d \n", prev_boundary);
                                fprintf(stderr, "\t****** >> lclg_node_idx: %d,  unispg_node_idx: %d\n", lclg_node_idx, unispg_node_idx);
                                fprintf(stderr, "\t****** >> lclg_node->start: %u,  lclg_node->end: %u\n", lclg_node->start, lclg_node->end);
                                fprintf(stderr, "\t****** >> unispg_node->start: %d,  unispg_node->end: %d\n", unispg_node->start, unispg_node->end);


                                CGraphnodeUnispg* node;

                                bool lclg_node_move = false;
                                bool unispg_node_move = false;


                                // If any of the node id is -1 => it is a fake node => go to boundary case check
                                // if (lclg_node->nodeid!=-1 && unispg_node->nodeid!=-1) {
                                    /******************************
                                    ** Creating graph nodes before boundary cases.
                                    *******************************/
                                    if (unispg_node->start < lclg_node->start) {
                                        // AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                        uint node_start_pos = 0;
                                        if (prev_boundary < unispg_node->start) {
                                            node_start_pos = unispg_node->start;
                                        } else if (prev_boundary == unispg_node->start) {
                                            node_start_pos = unispg_node->start;
                                        } else if (prev_boundary > unispg_node->start) {
                                            node_start_pos = prev_boundary;
                                        }
                                        if (unispg_node->end < lclg_node->start) {
                                            fprintf(stderr,"\t  ####  Graph node: ----------  |(s).................(e)|\n");
                                            // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            prev_boundary = unispg_node->end;
                                            MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                        } else if (unispg_node->end == lclg_node->start) {
                                            fprintf(stderr,"\t  ####  Graph node: ----------|(s).................(e)|\n");
                                            // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_S);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            prev_boundary = unispg_node->end;
                                            MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                        } else if (lclg_node->start < unispg_node->end) {
                                            // AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->start, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            if (unispg_node->end < lclg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: ------|(s)----.............(e)|\n");
                                                // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;

                                                prev_boundary = unispg_node->end;
                                                MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                            } else if (unispg_node->end == lclg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: ------|(s)----------(e)|\n");
                                                // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;


                                                prev_boundary = unispg_node->end;
                                                MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                                MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                            } else if (lclg_node->end < unispg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: -----|(s)----------(e)|-----\n");
                                                // AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;

                                                prev_boundary = lclg_node->end;
                                                MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                            }
                                        }
                                    } else if (unispg_node->start == lclg_node->start) {
                                        // AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S_LCLG_S);
                                        uint node_start_pos = 0;
                                        if (prev_boundary < unispg_node->start) {
                                            node_start_pos = unispg_node->start;
                                        } else if (prev_boundary == unispg_node->start) {
                                            node_start_pos = unispg_node->start;
                                        } else if (prev_boundary > unispg_node->start) {
                                            node_start_pos = prev_boundary;
                                        }


                                        if (unispg_node->end < lclg_node->end) {
                                            fprintf(stderr,"\t  ####  Graph node: |(s)----------.......(e)|\n");
                                            // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            prev_boundary = unispg_node->end;
                                            MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                        } else if (unispg_node->end == lclg_node->end) {
                                            fprintf(stderr,"\t  ####  Graph node: |(s)------------(e)|\n");
                                            // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            prev_boundary = unispg_node->end;
                                            MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                            MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                        } else if (lclg_node->end < unispg_node->end) {
                                            fprintf(stderr,"\t  ####  Graph node: |(s)---------(e)|---\n");
                                            // AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            prev_boundary = lclg_node->end;
                                            MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                        }
                                    
                                        
                                    } else if (lclg_node->start < unispg_node->start) {
                                        // AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
                                        uint node_start_pos = 0;
                                        if (prev_boundary < lclg_node->start) {
                                            node_start_pos = lclg_node->start;
                                        } else if (prev_boundary == lclg_node->start) {
                                            node_start_pos = lclg_node->start;
                                        } else if (prev_boundary > lclg_node->start) {
                                            node_start_pos = prev_boundary;
                                        }
                                        if (unispg_node->start < lclg_node->end) {
                                            // AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                            GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                            GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                            GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                            // Add new node
                                            node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->start, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                            new_no2gnode_unispg->Add(node);
                                            new_unispg_nodeid += 1;

                                            if (unispg_node->end < lclg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)....----------....(e)|\n");
                                                // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;

                                                prev_boundary = unispg_node->end;
                                                MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                            } else if (unispg_node->end == lclg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)....----------(e)|\n");
                                                // AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;

                                                prev_boundary = unispg_node->end;
                                                MoveUnispgNode(unispg_is_lastnode, unispg_node_move);
                                                MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                            } else if (lclg_node->end < unispg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)....----------(e)|---\n");
                                                // AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                                
                                                prev_boundary = lclg_node->end;
                                                MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                            }
                                        } else if (lclg_node->end == unispg_node->start) {
                                            fprintf(stderr,"\t  ####  Graph node: |(s)..........(e)|----------\n");
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
                                                node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                            }

                                            prev_boundary = lclg_node->end;
                                            MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                        } else if (lclg_node->end < unispg_node->start) {
                                            fprintf(stderr,"\t  ####  Graph node: |(s)..........(e)|  ----------\n");
                                            // AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
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
                                                node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                            }

                                            prev_boundary = lclg_node->end;
                                            MoveLclgNode(lclg_is_lastnode, lclg_node_move);
                                            // if (lclg_is_lastnode) {
                                            //     AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
                                            //     AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
                                            // }
                                        }
                                    }
                                // }




                                // if (lclg_is_lastnode && !unispg_is_lastnode) {
                                //     lclg_i+=1;
                                //     int lclg_node_idx = 1;
                                // } else if (!lclg_is_lastnode && unispg_is_lastnode) {
                                //     unispg_i+=1;
                                //     int unispg_node_idx = 1;
                                // } else if (!lclg_is_lastnode && !unispg_is_lastnode) {

                                // } else if (lclg_is_lastnode && unispg_is_lastnode) {
                                //     //  Tracking graph id.
                                //     lclg_i+=1;
                                //     unispg_i+=1;
                                //     //  Tracking node id.
                                //     int lclg_node_idx = 1;
                                //     int unispg_node_idx = 1;
                                // }





                                /******************************
                                ** Move node index in graphs.
                                *******************************/
                                if (lclg_node_move) {
                                    /******************************
                                    ** lclg needs to move and its not the last node of the graph.
                                    *******************************/
                                    fprintf(stderr, "\tlclg_node_idx Move ! %d\n", lclg_node_idx);
                                    lclg_node_idx += 1;
                                    unispg_is_lastnode = false;
                                }
                                if (unispg_node_move) {
                                    /******************************
                                    ** unispg needs to move and its not the last node of the graph.
                                    *******************************/
                                    fprintf(stderr, "\tunispg_node_idx Move ! %d\n", unispg_node_idx);
                                    unispg_node_idx += 1;
                                    lclg_is_lastnode = false;
                                }
                                if (!lclg_node_move && !unispg_node_move) {
                                    fprintf(stderr, "\t>>>> Local & global graph node cannot move\n");
                                    /******************************
                                    ** 1. Need to move to next node in lclg, but it's the last node of lclg
                                    ** 2. Need to move to next node in unispg, but it's the last node of unispg
                                    ** 3. Need to move to next node in both unispg & lclg, but both are last nodes.
                                    *******************************/
                                    /******************************
                                    ** When either lclg or unispg reaches its last node, 
                                    **  * 1. Document the current idx of the graph that does not reach the last node.
                                    **  * 2. For the graph that reaches the last node, move to the next graph.
                                    *******************************/
                                    if (lclg_is_lastnode && !unispg_is_lastnode) {
                                        fprintf(stderr, "\t\t>>>> lclg_is_lastnode && !unispg_is_lastnode\n");
                                        // 1. Need to move lclg, but it's the last node of lclg
                                        //     ####  Graph node: -----|(s)----------(e)|-----\n"
                                        //     ####  Graph node: |(s)---------(e)|---\n"
                                        //     ####  Graph node: |(s)....----------(e)|---\n"
                                        //     ####  Graph node: |(s)..........(e)|----------\n"
                                        //     ####  Graph node: |(s)..........(e)|  ----------\n"
                                        // uint node_start_pos = 0;
                                        if (prev_boundary < lclg_node->end) {
                                            // it is impossible. Just the insane check
                                            GError("Wrong boundaries. Check the code!!");
                                        } else if (prev_boundary == lclg_node->end) {
                                            if (unispg_node->start < prev_boundary) {
                                                // ####  Graph node: -----|(s)----------(e)|-----\n"
                                                // ####  Graph node: |(s)---------(e)|---\n"
                                                // ####  Graph node: |(s)....----------(e)|---\n"
                                                // prev_boundary remains the 'lclg_node->end'
                                            } else if (unispg_node->start >= prev_boundary) {
                                                // ####  Graph node: |(s)..........(e)|----------\n"
                                                // ####  Graph node: |(s)..........(e)|  ----------\n"
                                                prev_boundary = lclg_node->end;
                                            }
                                        } else if (prev_boundary > lclg_node->end) {
                                            // This situation has happend at least one. Insert the new unispg.
                                            prev_boundary = unispg_node->start;
                                        }
                                        // DO NOT CREATE NEW NODE. Compare the the current unispg node & boundary to the first node of the next lclg.
                                        // Move to the first node of the next lclg.
                                        lclg_next = true;
                                        // lclg_i+=1;
                                        // lclg_node_idx = 1;

                                        // node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                        // new_no2gnode_unispg->Add(node);
                                        // new_unispg_nodeid += 1;
                                    } else if (!lclg_is_lastnode && unispg_is_lastnode) {
                                        fprintf(stderr, "\t\t>>>> !lclg_is_lastnode && unispg_is_lastnode\n");
                                        //  2. Need to move unispg, but it's the last node of unispg
                                        //     ####  Graph node: ----------  |(s).................(e)|\n"
                                        //     ####  Graph node: ----------|(s).................(e)|\n"
                                        //     ####  Graph node: ------|(s)----.............(e)|\n"
                                        //     ####  Graph node: |(s)----------.......(e)|\n"
                                        //     ####  Graph node: |(s)....----------....(e)|\n"
                                        if (prev_boundary < unispg_node->end) {
                                            // it is impossible. Just the insane check
                                            GError("Wrong boundaries. Check the code!!");
                                        } else if (prev_boundary == unispg_node->end) {
                                            if (lclg_node->start < prev_boundary) {
                                                // ####  Graph node: ------|(s)----.............(e)|\n"
                                                // ####  Graph node: |(s)----------.......(e)|\n"
                                                // ####  Graph node: |(s)....----------....(e)|\n"
                                                // prev_boundary remains the 'unispg_node->end'
                                            } else if (lclg_node->start >= prev_boundary) {
                                                // ####  Graph node: ----------  |(s).................(e)|\n"
                                                // ####  Graph node: ----------|(s).................(e)|\n"
                                                prev_boundary = unispg_node->end;
                                            }
                                        } else if (prev_boundary > unispg_node->end) {
                                            // This situation has happend at least one. Insert the new unispg.
                                            prev_boundary = lclg_node->start;
                                        }
                                        // DO NOT CREATE NEW NODE. Compare the the current lclg node & boundary to the first node of the next unispg.
                                        // Move to the first node of the next unispg.
                                        unispg_next = true;
                                        // unispg_i+=1;
                                        // unispg_node_idx = 1;

                                        // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                        // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                        // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                        // // Add new node
                                        // fprintf(stderr, "^^^^ node_start_pos: %d,  unispg_node->end: %d \n", node_start_pos, lclg_node->end);
                                        // node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                        // new_no2gnode_unispg->Add(node);
                                        // new_unispg_nodeid += 1;
                                    } else if (!lclg_is_lastnode && !unispg_is_lastnode) {
                                        // This is an insane check. If any of unispg or lclg is not end, at least 1 must move. 
                                            GError("Wrong boundaries. Check the code!!");
                                    } else if (lclg_is_lastnode && unispg_is_lastnode) {
                                        fprintf(stderr, "\t\t>>>> lclg_is_lastnode && unispg_is_lastnode\n");
                                        // 3. Need to move both unispg & lclg, but both are last nodes.
                                        //     ####  Graph node: ------|(s)----------(e)|\n"
                                        //     ####  Graph node: |(s)------------(e)|\n"
                                        //     ####  Graph node: |(s)....----------(e)|\n"
                                        //  prev_boundary must be 'unispg_node->end && lclg_node->end'.
                                        if (unispg_node->start < lclg_node->start) {
                                            if (unispg_node->end < lclg_node->start) {
                                                fprintf(stderr,"\t  ####  Graph node: ----------  |(s).................(e)|\n");
                                                prev_boundary = unispg_node->end;
                                                unispg_next = true;
                                                // unispg_i+=1;
                                                // int unispg_node_idx = 1;

                                                // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                // // Add new node
                                                // node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                // new_no2gnode_unispg->Add(node);
                                                // new_unispg_nodeid += 1;
                                            } else if (unispg_node->end == lclg_node->start) {
                                                fprintf(stderr,"\t  ####  Graph node: ----------|(s).................(e)|\n");
                                                prev_boundary = unispg_node->end;
                                                unispg_next = true;
                                                // unispg_i+=1;
                                                // int unispg_node_idx = 1;

                                                // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                // // Add new node
                                                // node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                // new_no2gnode_unispg->Add(node);
                                                // new_unispg_nodeid += 1;
                                            } else if (lclg_node->start < unispg_node->end) {
                                                if (unispg_node->end < lclg_node->end) {
                                                    fprintf(stderr,"\t  ####  Graph node: ------|(s)----.............(e)|\n");
                                                    prev_boundary = unispg_node->end;
                                                    unispg_next = true;
                                                    // unispg_i+=1;
                                                    // int unispg_node_idx = 1;

                                                    // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                    // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                    // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                    // // Add new node
                                                    // node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                    // new_no2gnode_unispg->Add(node);
                                                    // new_unispg_nodeid += 1;
                                                } else if (unispg_node->end == lclg_node->end) {
                                                    fprintf(stderr,"\t  ####  Graph node: ------|(s)----------(e)|\n");
                                                    // Do not need to create a new node
                                                    prev_boundary = unispg_node->end;
                                                    unispg_next = true;
                                                    lclg_next = true;
                                                    // lclg_i+=1;
                                                    // unispg_i+=1;
                                                    // int lclg_node_idx = 1;
                                                    // int unispg_node_idx = 1;
                                                } else if (lclg_node->end < unispg_node->end) {
                                                    fprintf(stderr,"\t  ####  Graph node: -----|(s)----------(e)|-----\n");
                                                    prev_boundary = lclg_node->end;   
                                                    lclg_next = true;
                                                    // lclg_i+=1;
                                                    // int lclg_node_idx = 1;

                                                    // // Add new node
                                                    // node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                    // new_no2gnode_unispg->Add(node);
                                                    // new_unispg_nodeid += 1;
                                                }
                                            }
                                        } else if (unispg_node->start == lclg_node->start) {
                                            if (unispg_node->end < lclg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)----------.......(e)|\n");
                                                prev_boundary = unispg_node->end;   
                                                unispg_next = true;
                                                // unispg_i+=1;
                                                // int unispg_node_idx = 1;

                                                // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                // // Add new node
                                                // node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                // new_no2gnode_unispg->Add(node);
                                                // new_unispg_nodeid += 1;
                                            } else if (unispg_node->end == lclg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)------------(e)|\n");
                                                prev_boundary = unispg_node->end;
                                                unispg_next = true;
                                                lclg_next = true;
                                                // lclg_i+=1;
                                                // unispg_i+=1;
                                                // int lclg_node_idx = 1;
                                                // int unispg_node_idx = 1;
                                            } else if (lclg_node->end < unispg_node->end) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)---------(e)|---\n");
                                                prev_boundary = lclg_node->end;   
                                                lclg_next = true;
                                                // lclg_i+=1;
                                                // int lclg_node_idx = 1;

                                                // Add new node
                                                // node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                // new_no2gnode_unispg->Add(node);
                                                // new_unispg_nodeid += 1;
                                            }
                                        } else if (lclg_node->start < unispg_node->start) {
                                            if (unispg_node->start < lclg_node->end) {
                                                if (unispg_node->end < lclg_node->end) {
                                                    fprintf(stderr,"\t  ####  Graph node: |(s)....----------....(e)|\n");
                                                    prev_boundary = unispg_node->end;   
                                                    unispg_next = true;
                                                    // unispg_i+=1;
                                                    // int unispg_node_idx = 1;

                                                    // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                    // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                    // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                    // // Add new node
                                                    // node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                    // new_no2gnode_unispg->Add(node);
                                                    // new_unispg_nodeid += 1;
                                                } else if (unispg_node->end == lclg_node->end) {
                                                    fprintf(stderr,"\t  ####  Graph node: |(s)....----------(e)|\n");
                                                    prev_boundary = unispg_node->end;
                                                    unispg_next = true;
                                                    lclg_next = true;
                                                    // lclg_i+=1;
                                                    // unispg_i+=1;
                                                    // int lclg_node_idx = 1;
                                                    // int unispg_node_idx = 1;
                                                } else if (lclg_node->end < unispg_node->end) {
                                                    fprintf(stderr,"\t  ####  Graph node: |(s)....----------(e)|---\n");
                                                    prev_boundary = lclg_node->end;   
                                                    lclg_next = true;
                                                    // lclg_i+=1;
                                                    // int lclg_node_idx = 1;

                                                    // Add new node
                                                    // node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                    // new_no2gnode_unispg->Add(node);
                                                    // new_unispg_nodeid += 1;
                                                }
                                            } else if (lclg_node->end == unispg_node->start) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)..........(e)|----------\n");
                                                prev_boundary = lclg_node->end;   
                                                lclg_next = true;
                                                // lclg_i+=1;
                                                // int lclg_node_idx = 1;

                                                // Add new node
                                                // node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                // new_no2gnode_unispg->Add(node);
                                                // new_unispg_nodeid += 1;
                                            } else if (lclg_node->end < unispg_node->start) {
                                                fprintf(stderr,"\t  ####  Graph node: |(s)..........(e)|  ----------\n \n");
                                                prev_boundary = lclg_node->end;   
                                                lclg_next = true;
                                                // lclg_i+=1;
                                                // int lclg_node_idx = 1;

                                                // Add new node
                                                // node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                // new_no2gnode_unispg->Add(node);
                                                // new_unispg_nodeid += 1;
                                            }
                                        }
                                    }

                        // fprintf(stderr, "**** (%d) lclg_idx_start: %d, lclg_idx_end %d, lclg_limit: %d \n", lclg_i, lclg_idx_start, lclg_idx_end, lclg_limit);
                        // fprintf(stderr, "   >>> lclg boundary start: %d \n", lclg_start_pcs);
                        // fprintf(stderr, "   >>> lclg boundary end: %d \n", lclg_end_pcs);
                        // fprintf(stderr, "**** (%d) unispg_idx_start: %d, unispg_idx_end %d \n", unispg_i, unispg_idx_start, unispg_idx_end);
                        // fprintf(stderr, "   >>> unispg boundary start: %d \n", unispg_start_pcs);
                        // fprintf(stderr, "   >>> unispg boundary start: %d \n", unispg_end_pcs);

                                    /******************************
                                    ** Handling boundary cases.
                                    *******************************/
                                    if (lclg_i == lclg_idx_end && unispg_i < unispg_idx_end) {
                                        /******************************
                                        ** Move the lclg index. However, it has already been the last graph. (unispg_i should be unispg_idx_end-1)
                                        **  1. ####  Graph node: -----|(s)----------(e)|-----
                                        **  2. ####  Graph node: |(s)---------(e)|---
                                        **  3. ####  Graph node: |(s)....----------(e)|---
                                        **  4. ####  Graph node: |(s)..........(e)|----------
                                        **  5. ####  Graph node: |(s)..........(e)|  ----------
                                        *******************************/
                                        fprintf(stderr, "Move the lclg index. However, it has already been the last graph. (unispg_i should be unispg_idx_end-1)\n");
                                        if (lclg_node->nodeid==-1 && lclg_node->start != UINT_MAX && lclg_node->end != UINT_MAX) {
                                            if (lclg_node->end < unispg_node->start) {
                                                // Case 5
                                                node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                            } else {
                                                // Case 1, 2, 3, 4
                                                node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                            }
                                            prev_boundary = unispg_node->end;
                                        }
                                    } else if (lclg_i < lclg_idx_end && unispg_i == unispg_idx_end) {
                                        /******************************
                                        ** Move the unispg index. However, it has already been the last graph. (lclg_i should be lclg_idx_end-1)
                                        **   1. ####  Graph node: ----------  |(s).................(e)|
                                        **   2. ####  Graph node: ----------|(s).................(e)|
                                        **   3. ####  Graph node: ------|(s)----.............(e)|
                                        **   4. ####  Graph node: |(s)----------.......(e)|
                                        **   5. ####  Graph node: |(s)....----------....(e)|
                                        *******************************/
                                        fprintf(stderr, "Move the unispg index. However, it has already been the last graph. (lclg_i should be lclg_idx_end-1)\n");
                                        if (unispg_node->nodeid==-1 && unispg_node->start != UINT_MAX && unispg_node->end != UINT_MAX) {
                                            if (unispg_node->end < lclg_node->start) {
                                                // Case 1
                                                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                            } else {
                                                // Case 2, 3, 4, 5
                                                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                                                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                                                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                                                // Add new node
                                                node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov_s->Last(), lclg_node->capacity_s->Last(), 0);
                                                new_no2gnode_unispg->Add(node);
                                                new_unispg_nodeid += 1;
                                            }
                                            prev_boundary = lclg_node->end;
                                        }
                                    } else if (lclg_i == lclg_idx_end && unispg_i == unispg_idx_end) {
                                        /******************************
                                        ** Move the lclg & unispg indices. However, both of them have already reached the last graph.
                                        **  1. ####  Graph node: ------|(s)----------(e)|
                                        **  2. ####  Graph node: |(s)------------(e)|
                                        **  3. ####  Graph node: |(s)....----------(e)|
                                        *******************************/
                                        fprintf(stderr, "Move the lclg & unispg indices. However, both of them have already reached the last graph.\n");
                                    }
                                }
                            }
                        

                        }

                        if (unispg_next && unispg_i < unispg_idx_end) {
                            fprintf(stderr, "lclg_next: %d, lclg_i: %d, lclg_idx_start: %d, lclg_idx_end %d \n", lclg_next, lclg_i, lclg_idx_start, lclg_idx_end);
                            fprintf(stderr, "unispg_next: %d, unispg_i: %d, unispg_idx_start: %d, unispg_idx_end %d \n", unispg_next, unispg_i, unispg_idx_start, unispg_idx_end);
                            unispg_i += 1;
                            unispg_node_idx = 1;
                        }
                        if (lclg_next && lclg_i < lclg_idx_end) {
                            fprintf(stderr, "lclg_next: %d, lclg_i: %d, lclg_idx_start: %d, lclg_idx_end %d \n", lclg_next, lclg_i, lclg_idx_start, lclg_idx_end);
                            fprintf(stderr, "unispg_next: %d, unispg_i: %d, unispg_idx_start: %d, unispg_idx_end %d \n", unispg_next, unispg_i, unispg_idx_start, unispg_idx_end);
                            lclg_i += 1;
                            lclg_node_idx = 1;
                        }
                    }



                    GVec<bool>* is_passed_s_sink = new GVec<bool>(sample_num-1, false);
                    GVec<float>* cov_s_sink = new GVec<float>(sample_num-1, 0.0f);
                    GVec<float>* capacity_s_sink = new GVec<float>(sample_num-1, 0.0f);
                    CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, new_unispg_nodeid, is_passed_s_sink, cov_s_sink, capacity_s_sink, false, 0, 0, 0);
                    new_no2gnode_unispg->Add(sink);


                    // Printing lclg & unispg that are gonna be processed.
                    // fprintf(stderr, "Inside `process_ovp_graphs`!!!\n");
                    // for (int unispg_i = unispg_idx_start; unispg_i < unispg_idx_end; unispg_i++) {
                    //     fprintf(stderr, "unispg_index -> unispg_i: %d\n", unispg_i);
                    //     for (int unispg_j = 0; unispg_j < no2gnode_unispg[s][unispg_i].Count(); unispg_j++) {
                    //         fprintf(stderr, "** no2gnode_unispg[i][j]->start: %d\n", no2gnode_unispg[s][unispg_i][unispg_j]->start);
                    //         fprintf(stderr, "** no2gnode_unispg[i][j]->end: %d\n", no2gnode_unispg[s][unispg_i][unispg_j]->end);
                    //     }
                    // }
                    // for (int lclg_i = lclg_idx_start; lclg_i < lclg_idx_end; lclg_i++) {
                    //     fprintf(stderr, "lclg_index -> lclg_i: %d\n", lclg_i);
                    //     for (int lclg_j = 0; lclg_j < no2gnode[lclg_i].Count(); lclg_j++) {
                    //         fprintf(stderr, "** no2gnode[i][j]->start: %d\n", no2gnode[lclg_i][lclg_j]->start);
                    //         fprintf(stderr, "** no2gnode[i][j]->end: %d\n", no2gnode[lclg_i][lclg_j]->end);
                    //     }
                    // }

    // GStr node_g(g);
    GStr strand_symbol;
    if (s == 0) {
        strand_symbol = "-";
    } else if (s == 1) {
        strand_symbol = "+";
    }
    // GStr bundle_start("");
    // GStr bundle_end("");
    /****************
     **  Writing out the visualization graph for the global graph.
     ****************/

    // fprintf(stderr, "&& unispg_gp->current_gidx: %d\n", unispg_gp->current_gidx[s]-1);
    // GPVec<CGraphnode>** unispg_gp->no2gnode_unispg = unispg_gp->get_no2gnodeGp();
    fprintf(stderr, "unispg_gp->no2gnode_unispg[s][unispg_gp->current_gidx[s]-1].Count(): %d \n", new_no2gnode_unispg->Count());
    if (new_no2gnode_unispg->Count() > 1) {
        // int refstart_unispg = unispg_gp->no2gnode_unispg[s][g][1]->start;
        // int refend_unispg = unispg_gp->no2gnode_unispg[s][g][unispg_gp->no2gnode_unispg[s][g].Count()-2]->end;
        // GVec<int>* graphno_unispg = unispg_gp->get_graphnoGp();
        // GVec<int>* edgeno_unispg = unispg_gp->get_edgenoGp();

        fprintf(stderr,"Traversing the universal splice graph!!!\n");
        // fprintf(stderr,"Unispg %d_%d_%d_%d {", bdata->start, bdata->end, s, g);
        // graphno[s][b]: number of nodes in graph.
        fprintf(stderr,"new_no2gnode_unispg->Count(): %d !!!\n", new_no2gnode_unispg->Count());
        for(int nd=0;nd<new_no2gnode_unispg->Count();nd++) {
            fprintf(stderr,"%d[start=%d end=%d];",nd,new_no2gnode_unispg->Get(nd)->start,new_no2gnode_unispg->Get(nd)->end);
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
            // GStr node_name = "Node_" + bundle_start + "_" + bundle_end + "_" + node_g + "_" + node_nd;
            GStr node_name = "Node_" + node_nd;
            fprintf(stderr, "node_name: %s\n", node_name.chars());

            if (nd == 0) {
                node_start = new_no2gnode_unispg->Get(1)->start-20;
                node_end = new_no2gnode_unispg->Get(1)->start;
            } else if (nd == new_no2gnode_unispg->Count()-1){
                node_start = new_no2gnode_unispg->Get(new_no2gnode_unispg->Count()-2)->end-1;
                node_end = 	new_no2gnode_unispg->Get(new_no2gnode_unispg->Count()-2)->end+20;
            } else {
                node_start = new_no2gnode_unispg->Get(nd)->start-1;
                node_end = new_no2gnode_unispg->Get(nd)->end;		
            }


            if (nd == 0) {
                if(s == 0) {
                    fprintf(node_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                } else if (s == 1) {
                    fprintf(node_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                }
            } else if (nd == new_no2gnode_unispg->Count()-1){
                if(s == 0) {
                // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\tNODE\t%f\t+\n", lclg_nonoverlap[s][g][lclg_nonoverlap[s][g].Count()-2]->end, lclg_nonoverlap[s][g][lclg_nonoverlap[s][g].Count()-2]->end+200, 0);

                    fprintf(node_cov_neg_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t+\n", node_start, node_end, node_name.chars());
                } else if (s == 1) {
                    // fprintf(node_cov_pos_bed, "chr22\t%d\t%d\tNODE\t%f\t-\n", lclg_nonoverlap[s][g][lclg_nonoverlap[s][g].Count()-2]->end, lclg_nonoverlap[s][g][lclg_nonoverlap[s][g].Count()-2]->end+200, 0);

                    fprintf(node_cov_pos_bed_unispg_vec.Get(fidx), "chr22\t%d\t%d\t%s\t0\t-\n", node_start, node_end, node_name.chars());
                }
            } else {
                if(s == 0) {
                    // fprintf(node_cov_neg_bed, "chr22\t%d\t%d\t%f\t%s\n", lclg_nonoverlap[s][g][nd]->start, lclg_nonoverlap[s][g][nd]->end, lclg_nonoverlap[s][g][nd]->cov, strand_symbol.chars());
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
    }


                    
                    fprintf(stderr, "\n\n\n");

                    unispg_idx_start = current_gidx[s];
                    unispg_idx_end = current_gidx[s];

                    // lclg_idx = lclg_idx_end;
                    lclg_idx_start = lclg_idx;
                    lclg_idx_end = lclg_idx;
                }

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

void UnispgGp::MoveUnispgNode(bool& unispg_is_lastnode, bool& unispg_node_move) {
    if (!unispg_is_lastnode) {
        unispg_node_move = true;
    }
}

void UnispgGp::MoveLclgNode(bool& lclg_is_lastnode, bool& lclg_node_move) {
    if (!lclg_is_lastnode) {
        lclg_node_move = true;
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