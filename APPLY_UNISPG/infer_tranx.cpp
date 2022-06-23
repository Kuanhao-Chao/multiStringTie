#include "infer_tranx.h"

void infer_transcripts_unispg(BundleData* bundle, GPVec<UnispgGp>** graphs_vec) {

    count_good_junctions(bundle);
	
    // geneno = build_graphs_unispg(bundle, unispg);
    fprintf(stderr, "Inside `infer_transcripts_unispg`\n");
    fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
    fprintf(stderr, "bundle->bpcov[1].Count(): %d\n", bundle->bpcov[1].Count());
    fprintf(stderr, "bundle->start - end: %d - %d\n", bundle->start, bundle->end);


    for(int s=0;s<2;s+=1) {
        fprintf(stderr, "\t@@ graphs_vec[%d]: %d\n", s, graphs_vec[s]->Count());
        for (int i=0; i<graphs_vec[s]->Count(); i++) {
            fprintf(stderr, "\t\t graphs_vec bound: %u - %u \n", graphs_vec[s]->Get(i)->get_refstart(), graphs_vec[s]->Get(i)->get_refend());
        }
    }


    uint unispg_nstart = 0;
    uint unispg_nend = 0;

    int refstart = bundle->start;
    int refend = bundle->end;
    float node_coverage = 0.0;
    float node_coverage_neg = 0.0;
    float node_coverage_uns = 0.0;
    float node_coverage_pos = 0.0;

	/*****************************
	 ** Calculate bundle pos / neg coverage ratio
	 *****************************/
    float bundle_coverage_neg = get_cov(0, 0, refend-refstart+1, bundle->bpcov);
    float bundle_coverage_uns = get_cov(1, 0, refend-refstart+1, bundle->bpcov);
    float bundle_coverage_pos = get_cov(2, 0, refend-refstart+1, bundle->bpcov);

    float bundle_coverage_ratio_neg = get_cov(0, 0, refend-refstart+1, bundle->bpcov);
    float bundle_coverage_ratio_pos = get_cov(0, 0, refend-refstart+1, bundle->bpcov);

    if (bundle_coverage_neg == 0.0 && bundle_coverage_pos == 0.0) {
        bundle_coverage_ratio_neg = 0.5;
        bundle_coverage_ratio_pos = 0.5;
    } else {
        bundle_coverage_ratio_neg = bundle_coverage_neg/(bundle_coverage_neg+bundle_coverage_pos);
        bundle_coverage_ratio_pos = 1-bundle_coverage_ratio_neg;
    }

    fprintf(stderr, "bundle_coverage_neg: %f (ratio: %f); bundle_coverage_uns: %f; bundle_coverage_pos: %f (ratio: %f)\n", bundle_coverage_neg, bundle_coverage_ratio_neg, bundle_coverage_uns, bundle_coverage_pos, bundle_coverage_ratio_pos);

    if (graphs_vec[0]->Count() == 0) {
        fprintf(stderr, "####### Positive strand!!!\n");
        for (int g=0; g<graphs_vec[1]->Count(); g++) {
            for (int n=1; n<graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Count()-1; n++) {
                
                fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[1]->Get(g)->get_refstart(), graphs_vec[1]->Get(g)->get_refend(), g, n, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->nodeid, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end);

                unispg_nstart = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start;
                unispg_nend = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end;
                
                // node_coverage = get_cov(2, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

                node_coverage = get_cov_sign(2, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

                fprintf(stderr, "g: %d; n: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\n", g, n, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->nodeid, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end, node_coverage, node_coverage/(unispg_nend-unispg_nstart+1), graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s->Get(0), graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s->Get(0)/(unispg_nend-unispg_nstart+1));
                fprintf(stderr, "node start - end: %u - %u\n", graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end);

            }
        }
    } else if (graphs_vec[1]->Count() == 0) {
        fprintf(stderr, "####### Negative strand!!!\n");
        for (int g=0; g<graphs_vec[0]->Count(); g++) {
            for (int n=1; n<graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Count()-1; n++) {
                
                fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[0]->Get(g)->get_refstart(), graphs_vec[0]->Get(g)->get_refend(), g, n, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->nodeid, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end);

                unispg_nstart = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start;
                unispg_nend = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end;
                
                // node_coverage = get_cov(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

                node_coverage = get_cov_sign(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);
                fprintf(stderr, "g: %d; n: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\n", g, n, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->nodeid, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end, node_coverage, node_coverage/(unispg_nend-unispg_nstart+1), graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s->Get(0), graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s->Get(0)/(unispg_nend-unispg_nstart+1));
            }
        }
    } else {
        // Positive & negative graphs overlap!!
        int pos_g_idx = 0;
        int neg_g_idx = 0;
        int pos_n_idx = 1;
        int neg_n_idx = 1;
        
        fprintf(stderr, "graphs_vec[1]->Count(): %d;  graphs_vec[0]->Count(): %d\n", graphs_vec[1]->Count(), graphs_vec[0]->Count());

        while (pos_g_idx < graphs_vec[1]->Count() || neg_g_idx < graphs_vec[0]->Count()) {
            if ((pos_g_idx >= graphs_vec[1]->Count()) || (neg_g_idx >= graphs_vec[0]->Count())) {

                break;
                if ((pos_g_idx >= graphs_vec[1]->Count()) && (neg_g_idx >= graphs_vec[0]->Count())) {
                    break;
                } else {
                    if (pos_g_idx >= graphs_vec[1]->Count()) {
                        pos_g_idx = graphs_vec[1]->Count() - 1;
                        // neg_g_idx++;
                        // neg_n_idx = 1;
                    }
                    if (neg_g_idx >= graphs_vec[0]->Count()) {
                        neg_g_idx = graphs_vec[0]->Count()-1;
                        // pos_g_idx++;
                        // pos_n_idx = 1;
                    }
                }
            }
            fprintf(stderr, "pos_g_idx: %d;  neg_g_idx: %d\n", pos_g_idx, neg_g_idx);

            int pos_g_start = graphs_vec[1]->Get(pos_g_idx)->refstart;
            int pos_g_end = graphs_vec[1]->Get(pos_g_idx)->refend;

            int neg_g_start = graphs_vec[0]->Get(neg_g_idx)->refstart;
            int neg_g_end = graphs_vec[0]->Get(neg_g_idx)->refend;


            bool g_overlap = segs_overlap(pos_g_start, pos_g_end, neg_g_start, neg_g_end);
            if (g_overlap) {
                // Process the overlapping graphs

                fprintf(stderr, "graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1: %d;  graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1: %d\n", graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count(), graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count());
                while ((pos_n_idx < graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1) || (neg_n_idx < graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1)) {

                    if ((pos_n_idx == graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1) || (neg_n_idx == graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1)) {
                        if (pos_n_idx == graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1) {
                            // Positive graph reach the last node.
                            pos_g_idx++;
                            pos_n_idx = 1;
                        }
                        if (neg_n_idx == graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1) {
                            neg_g_idx++;
                            neg_n_idx = 1;
                        }
                        break;
                    }

                    uint pos_n_start = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->start;
                    uint pos_n_end = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->end;
                    uint neg_n_start = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->start;
                    uint neg_n_end = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->end;
                    bool n_overlap = segs_overlap(pos_n_start, pos_n_end, neg_n_start, neg_n_end);
                    if (n_overlap) {
                        // Process the overlapping nodes


                        // uint overlap_n_start = pos_n_start<neg_n_start ? pos_n_start : neg_n_start;
                        // uint overlap_n_end = pos_n_end>neg_n_end ? pos_n_end : neg_n_end;
                        // This is the the unstranded reads coverage
                        // node_coverage_neg = get_cov(0, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
                        node_coverage_uns = get_cov(1, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
                        // node_coverage_pos = get_cov(2, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
                        float node_coverage_pos = 0.0;
                        float node_coverage_neg = 0.0;

                        calculate_ovp_coverage(pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend, bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, node_coverage_pos, node_coverage_neg, bundle->bpcov);

                        // node_coverage = get_cov(1, overlap_n_start-refstart, overlap_n_end-refstart, bundle->bpcov);

                        fprintf(stderr, "pos_g_idx: %d; pos_n_idx: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\nneg_g_idx: %d; neg_n_idx: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\nunstranded node_coverage (%f)\n", pos_g_idx, pos_n_idx, graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->nodeid, pos_n_start, pos_n_end, node_coverage_pos, node_coverage_pos/(pos_n_end-pos_n_start), graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->cov_s->Get(0), graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->cov_s->Get(0)/(pos_n_end-pos_n_start), neg_g_idx, neg_n_idx, graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->nodeid, neg_n_start, neg_n_end, node_coverage_neg, node_coverage_neg/(neg_n_end-neg_n_start), graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->cov_s->Get(0), graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->cov_s->Get(0)/(neg_n_end-neg_n_start), node_coverage_uns);

                        if (pos_n_end < neg_n_end) {
                            // positive ends first
                            pos_n_idx++;
                        } else if (neg_n_end < pos_n_end) {
                            // negative ends first
                            neg_n_idx++;
                        } else if (neg_n_end == pos_n_end) {
                            pos_n_idx++;
                            neg_n_idx++;
                        }
                    } else {
                        if (pos_n_end < neg_n_start) {
                            // Positive is in the front & not overlap
                            pos_n_idx++;               
                        } else if (neg_n_end < pos_n_start) {
                            // negative is in the front & not overlap
                            neg_n_idx++;
                        }
                    }
                }

                if (pos_g_end < neg_g_end) {
                    // positive ends first
                    pos_g_idx++;
                } else if (neg_g_end < pos_g_end) {
                    // negative ends first
                    neg_g_idx++;
                } else if (neg_g_end == pos_g_end) {
                    pos_g_idx++;
                    neg_g_idx++;
                }
            } else {
                if (pos_g_end < neg_g_start) {
                    // Positive is in the front & not overlap
                    pos_g_idx++;               
                } else if (neg_g_end < pos_g_start) {
                    // negative is in the front & not overlap
                    neg_g_idx++;
                }
            }


            

        }
    }
    // graphs_vec[0]->Count();
    // graphs_vec[1]->Count();


    // Now I need to determine the overlapping condition of positive & negative graphs.

    // for (int i=0; i<unispg->no2gnode_unispg[1][0].Count(); i++) {
    //     fprintf(stderr, "unispg[%d] nodeid: %d\n", i, i<unispg->no2gnode_unispg[1][0].Get(i)->nodeid);
    // }

	// /*****************************
	//  ** Applying reads to universal graphs.
	//  *****************************/
    // GList<CReadAln>& readlist = bundle->readlist;
    // // GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
    // // Create here!!
    // GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    // // Create here!!
    // CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
    // // Create here!!
    // GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    // // Create here!!
    // GVec<int> lastgpos[2];
    // // Input variable!!
    // GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    // // Input variable!!
    // GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    // // GVec<int> trnumber[2]; // how many transfrags are on a strand s, in a graph g -> I can find out this from transfrag[s][g].Count()
    // // int ngraph[2]={0,0};   // how many graphs are in each strand: negative (0), or positive(1) -> keep one for each bundle


    // // Input variable!!
    // GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
    

    // int bno[2]={0,0};
    
    
    // // for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
    //     // int s=sno/2; // adjusted strand due to ignoring neutral strand
    // int s = 1;
    // //char strnd='-';
    // //if(s) strnd='+';

    // // bundle2graph[s]=NULL;
    // // if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    // transfrag[s]=NULL;
    // no2gnode[s]=NULL;
    // tr2no[s]=NULL;
    // gpos[s]=NULL;

    // uint unispg_nstart = 0;
    // uint unispg_nend = 0;

    // int refstart = unispg ->get_refstart();
    // int refend = unispg->get_refend();

    // int readlist_idx = 0;
    // float node_coverage = 0.0;
    // fprintf(stderr, "unispg->graphno_unispg[s]: %d\n", unispg->graphno_unispg[s]);
    // fprintf(stderr, "unispg->edgeno_unispg[s]: %d\n", unispg->edgeno_unispg[s]);
    // for (int i=1; i<unispg->graphno_unispg[s]-1; i++) {

    //     fprintf(stderr, "%d; ref (%d - %d); unispg[%d] nodeid: %d (%u - %u)\n", i, refstart, refend, i, unispg->no2gnode_unispg[1][0].Get(i)->nodeid, unispg->no2gnode_unispg[1][0].Get(i)->start, unispg->no2gnode_unispg[1][0].Get(i)->end);




    //     unispg_nstart = unispg->no2gnode_unispg[1][0].Get(i)->start;
    //     unispg_nend = unispg->no2gnode_unispg[1][0].Get(i)->end;

    //     // readlist[readlist_idx]->start;
    //     // readlist[readlist_idx]->end;
        
    //     node_coverage = get_cov(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

    //     // node_coverage = get_cov_sign(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);
    //     fprintf(stderr, "unispg[%d] nodeid: %d (%u - %u) => node_coverage (%f / %f)\n", i, unispg->no2gnode_unispg[1][0].Get(i)->nodeid, unispg->no2gnode_unispg[1][0].Get(i)->start, unispg->no2gnode_unispg[1][0].Get(i)->end, node_coverage, unispg->no2gnode_unispg[1][0].Get(i)->cov_s->Get(0));
    // }
}
