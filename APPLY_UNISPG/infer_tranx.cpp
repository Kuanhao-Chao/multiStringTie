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

    /*****************************
     **    'get_fragment_pattern' function
     **        because of this going throu
     *****************************/
    for (int n=0; n < bundle->readlist.Count(); n++) {
    	float single_count=bundle->readlist[n]->read_count;
        // fprintf(stderr, ">> single_count: %f;  epsilon: %f\n", single_count, epsilon);
        for (int i=0; i<bundle->readlist[n]->pair_idx.Count(); i++) {
            // fprintf(stderr, ">> single_count: %d\n", i);
            int np = bundle->readlist[n]->pair_idx[i];
            if (np > -1) {
                single_count -= bundle->readlist[n]->pair_idx[i];
                fprintf(stderr, ">> n: %d;   np: %d\n", n, np);
                if (n < np) {
                    // fprintf(stderr, ">> n < np: %d\n", n < np);
    				get_fragment_pattern(bundle->readlist, n, np, graphs_vec);
                        // readlist,n,np,readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
                }
            }
        }
        if (single_count > epsilon) {
    		get_fragment_pattern(bundle->readlist, n, -1, graphs_vec);
                // readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
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
     ** Redistribute by whole bundleData
     *****************************/
    float bundle_coverage_neg = get_cov(0, 0, refend-refstart+1, bundle->bpcov);
    float bundle_coverage_uns = get_cov(1, 0, refend-refstart+1, bundle->bpcov);
    float bundle_coverage_pos = get_cov(2, 0, refend-refstart+1, bundle->bpcov);

    float bundle_coverage_ratio_neg = 0.0;
    // get_cov(0, 0, refend-refstart+1, bundle->bpcov);
    float bundle_coverage_ratio_pos = 0.0;
    // get_cov(2, 0, refend-refstart+1, bundle->bpcov);

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
                float target_cov_ratio = 0.0;
                fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[1]->Get(g)->get_refstart(), graphs_vec[1]->Get(g)->get_refend(), g, n, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->nodeid, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end);

                unispg_nstart = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start;
                unispg_nend = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end;
                
                // node_coverage = get_cov(2, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

                node_coverage = get_cov_sign(2, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

                fprintf(stderr, "g: %d; n: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\n", g, n, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->nodeid, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end, node_coverage, node_coverage/(unispg_nend-unispg_nstart+1), graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s->Get(0), graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s->Get(0)/(unispg_nend-unispg_nstart+1));
                fprintf(stderr, "node start - end: %u - %u\n", graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end);
                target_cov_ratio = node_coverage/(graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s->Get(0));
                ratio_file << target_cov_ratio << "\n";
                fprintf(stderr, "Single node!!\n");
            }


        }
    } else if (graphs_vec[1]->Count() == 0) {
        fprintf(stderr, "####### Negative strand!!!\n");
        for (int g=0; g<graphs_vec[0]->Count(); g++) {
            for (int n=1; n<graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Count()-1; n++) {
                float target_cov_ratio = 0.0;
                fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[0]->Get(g)->get_refstart(), graphs_vec[0]->Get(g)->get_refend(), g, n, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->nodeid, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end);

                unispg_nstart = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start;
                unispg_nend = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end;
                
                // node_coverage = get_cov(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

                node_coverage = get_cov_sign(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);
                fprintf(stderr, "g: %d; n: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\n", g, n, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->nodeid, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end, node_coverage, node_coverage/(unispg_nend-unispg_nstart+1), graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s->Get(0), graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s->Get(0)/(unispg_nend-unispg_nstart+1));
                target_cov_ratio = node_coverage/(graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s->Get(0));
                ratio_file << target_cov_ratio << "\n";
                fprintf(stderr, "Single node!!\n");
            }
        }
    } else {
        // Positive & negative graphs overlap!!
        int pos_g_idx = 0;
        int neg_g_idx = 0;
        int pos_n_idx = 1;
        int neg_n_idx = 1;

        int pos_g_num = graphs_vec[1]->Count();
        int neg_g_num = graphs_vec[0]->Count();
        int pos_n_num = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1;
        int neg_n_num = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1;

        bool pos_reach_end = false;
        bool neg_reach_end = false;
        bool pos_reach_end_first_time = true;
        bool neg_reach_end_first_time = true;

        int chaining_num = 1;

        int last_ovp_end = 0;
        char chaining_hold_strand = '.';
        float chaining_hold_cov = 0.0;

        // GPVec<CGraphnodeUnispg>* ovp_nodes = new GPVec<CGraphnodeUnispg>[40];


        // First pass => check if nodes are close enough to share pos / neg ratio
        // int fp_pos_g_idx = 0;
        // int fp_neg_g_idx = 0;
        // int fp_pos_n_idx = 1;
        // int fp_neg_n_idx = 1;


        // Second pass => actually process  nodes.
        while(!pos_reach_end || !neg_reach_end) {
            fprintf(stderr, "*** Current position:!!!\n");
            fprintf(stderr, "*** Positive g: %d,  n: %d\n", pos_g_idx, pos_n_idx);
            fprintf(stderr, "*** Negative g: %d,  n: %d\n", neg_g_idx, neg_n_idx);

            pos_n_num = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1;
            neg_n_num = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1;

            /*****************************
             ** Check if pos / neg reach end
             *****************************/
            if (pos_g_idx == pos_g_num-1 && pos_n_idx == pos_n_num-1) {
                pos_reach_end = true;
                if (pos_reach_end_first_time) {
                    // Get the coverage of the last positive node!!!
                    pos_reach_end_first_time = false;
                }
                // Get the cov of the last positive node
            }
            if (neg_g_idx == neg_g_num-1 && neg_n_idx == neg_n_num-1) {
                neg_reach_end = true;
                if (neg_reach_end_first_time) {
                    // Get the coverage of the last negative node!!!
                    neg_reach_end_first_time = false;
                }
                // fprintf(stderr, "Reaching the end of the negative graph!!!!\n");
                // Get the cov of the last negative node
            }

            if (pos_reach_end && neg_reach_end) break; 
            /*****************************
             ** Getting node position
            *****************************/
            uint pos_n_start = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->start;
            uint pos_n_end = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->end;
            uint neg_n_start = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->start;
            uint neg_n_end = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->end;

            // if (last_ovp_end == 0) {
            //     last_ovp_end = pos_n_start<neg_n_start ? pos_n_start : neg_n_start;
            // }


            // /*****************************
            //  ** Redistribute by 2 nodes coverage ratio
            //  *****************************/
            // uint lcl_boundary_start = min(pos_n_start, neg_n_start);
            // uint lcl_boundary_end = max(pos_n_end, neg_n_end);
            // // fprintf(stderr, "lcl_boundary_start: %u; lcl_boundary_end: %u \n", lcl_boundary_start, lcl_boundary_end);

            // float lcl_coverage_neg = get_cov(0, lcl_boundary_start-refstart, lcl_boundary_end-refstart+1, bundle->bpcov);
            // float lcl_coverage_uns = get_cov(1, lcl_boundary_start-refstart, lcl_boundary_end-refstart+1, bundle->bpcov);
            // float lcl_coverage_pos = get_cov(2, lcl_boundary_start-refstart, lcl_boundary_end-refstart+1, bundle->bpcov);

            // // fprintf(stderr, "lcl_coverage_neg: %f; bundle_coverage_uns: %f; lcl_coverage_pos: %f\n", lcl_coverage_neg, bundle_coverage_uns, lcl_coverage_pos);


            // bundle_coverage_ratio_neg = 0.0;
            // // get_cov(0, 0, refend-refstart+1, bundle->bpcov);
            // bundle_coverage_ratio_pos = 0.0;
            // // get_cov(2, 0, refend-refstart+1, bundle->bpcov);

            // if (lcl_coverage_neg == 0.0 && lcl_coverage_pos == 0.0) {
            //     bundle_coverage_ratio_neg = 0.5;
            //     bundle_coverage_ratio_pos = 0.5;
            // } else {
            //     bundle_coverage_ratio_neg = lcl_coverage_neg/(lcl_coverage_neg+lcl_coverage_pos);
            //     bundle_coverage_ratio_pos = 1-bundle_coverage_ratio_neg;
            // }

            // fprintf(stderr, "lcl_coverage_neg: %f (ratio: %f); bundle_coverage_uns: %f; lcl_coverage_pos: %f (ratio: %f)\n", lcl_coverage_neg, bundle_coverage_ratio_neg, bundle_coverage_uns, lcl_coverage_pos, bundle_coverage_ratio_pos);



            /*****************************
             ** Redistribute by 2 graphs coverage ratio
             *****************************/
            // uint pos_g_start = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(1)->start;
            // uint pos_g_end = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_num-1)->end;
            // uint neg_g_start = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(1)->start;
            // uint neg_g_end = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_num-1)->end;

            // uint lcl_boundary_start = min(pos_g_start, neg_g_start);
            // uint lcl_boundary_end = max(pos_g_end, neg_g_end);

            // float lcl_coverage_neg = get_cov(0, lcl_boundary_start-refstart, lcl_boundary_end-refstart+1, bundle->bpcov);
            // float lcl_coverage_uns = get_cov(1, lcl_boundary_start-refstart, lcl_boundary_end-refstart+1, bundle->bpcov);
            // float lcl_coverage_pos = get_cov(2, lcl_boundary_start-refstart, lcl_boundary_end-refstart+1, bundle->bpcov);

            // // fprintf(stderr, "lcl_coverage_neg: %f; bundle_coverage_uns: %f; lcl_coverage_pos: %f\n", lcl_coverage_neg, bundle_coverage_uns, lcl_coverage_pos);


            // bundle_coverage_ratio_neg = 0.0;
            // // get_cov(0, 0, refend-refstart+1, bundle->bpcov);
            // bundle_coverage_ratio_pos = 0.0;
            // // get_cov(2, 0, refend-refstart+1, bundle->bpcov);

            // if (lcl_coverage_neg == 0.0 && lcl_coverage_pos == 0.0) {
            //     bundle_coverage_ratio_neg = 0.5;
            //     bundle_coverage_ratio_pos = 0.5;
            // } else {
            //     bundle_coverage_ratio_neg = lcl_coverage_neg/(lcl_coverage_neg+lcl_coverage_pos);
            //     bundle_coverage_ratio_pos = 1-bundle_coverage_ratio_neg;
            // }

            // fprintf(stderr, "lcl_coverage_neg: %f (ratio: %f); bundle_coverage_uns: %f; lcl_coverage_pos: %f (ratio: %f)\n", lcl_coverage_neg, bundle_coverage_ratio_neg, bundle_coverage_uns, lcl_coverage_pos, bundle_coverage_ratio_pos);




            // bool n_overlap = segs_overlap_chain(pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend, bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, bundle->bpcov);

            bool n_overlap = segs_overlap(pos_n_start, pos_n_end, neg_n_start, neg_n_end);
            // bool n_overlap = segs_overlap(pos_n_start, pos_n_end, neg_n_start, neg_n_end);

            if (n_overlap) {
                /*****************************
                 ** Add the node that's in the back into the chain
                 *****************************/
                chaining_num += 1;
                fprintf(stderr, "\t@@@@@@@ Overlapping!!!! Chainning!!!!!!\n");
                fprintf(stderr, "\tCurrent chaining_num: %d\n", chaining_num);

                // if (pos_n_end < neg_n_end) {
                //     // ------......(e)|: I can get the coverage of the positive node.!!
                //     // CGraphnodeUnispg* ovp_node_neg = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx);
                //     // ovp_nodes->Add(ovp_node_neg);

                // } else if (pos_n_end > neg_n_end) {
                //     // ------(e)|-----: I can get the coverage of the negative node.!!
                //     // CGraphnodeUnispg* ovp_node_pos = graphs_vec[1]->Get(neg_g_idx)->no2gnode_unispg[1][0].Get(neg_n_idx);
                //     // ovp_nodes->Add(ovp_node_pos);

                // } else if (pos_n_end == neg_n_end) {
                //     // ------(e)|: I can get the coverage of both nodes.!!
                //     // CGraphnodeUnispg* ovp_node_pos = graphs_vec[1]->Get(neg_g_idx)->no2gnode_unispg[1][0].Get(neg_n_idx);
                //     // ovp_nodes->Add(ovp_node_pos);
                // }
            } else {
                /*****************************
                 ** start processing the chain
                 *****************************/
                fprintf(stderr, "&&&&&&& Stop Overlapping!!!! Stop chainning!!!!!\n");
                chaining_num = 1;

                /*****************************
                 ** empty the chain
                 *****************************/
                // ovp_nodes->Clear();
                
                /*****************************
                 ** Add the node that's in the front
                 *****************************/

                // if (pos_n_end < neg_n_start) {
                //     // ------  |(s)......: add the positive node
                //     // CGraphnodeUnispg* ovp_node_pos = graphs_vec[1]->Get(neg_g_idx)->no2gnode_unispg[1][0].Get(neg_n_idx);
                // } else if (pos_n_start > neg_n_end) {
                //     // .....(e)|  --------: add the negative node
                //     // CGraphnodeUnispg* ovp_node_neg = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx);
                // } else {
                //     // This is an error => cannot enter here.
                // }
                // CGraphnodeUnispg* ovp_node_pos = graphs_vec[1]->Get(neg_g_idx)->no2gnode_unispg[1][0].Get(neg_n_idx);
                // CGraphnodeUnispg* ovp_node_neg = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx);
            }


            float expected_cov_pos = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->cov_s->Get(0);
            float expected_cov_neg = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->cov_s->Get(0);

            fprintf(stderr, "\t>>> expected_cov_pos: %f;  expected_cov_neg: %f \n", expected_cov_pos, expected_cov_neg);

            float end_chaining_cov = 0.0;
            float target_cov_ratio = 0.0;
            /*****************************
             ** Push to the next node.
             *****************************/
            // I need to track which strand do I push back!!!
            // Calculate the coverage of the last node in the strand that's being pushed back.

            // fprintf(stderr, "neg_g_idx: %d; neg_g_num: %d; neg_n_idx: %d; neg_n_num: %d; neg_reach_end: %d\n", neg_g_idx, neg_g_num, neg_n_idx, neg_n_num, neg_reach_end);

            // fprintf(stderr, "pos_g_idx: %d; pos_g_num: %d; pos_n_idx: %d; pos_n_num: %d; pos_reach_end: %d\n", pos_g_idx, pos_g_num, pos_n_idx, pos_n_num, pos_reach_end);

            if (pos_n_end < neg_n_end) {
                // Push positive node
                if (!pos_reach_end) {
                    fprintf(stderr, "\t>>> Pushing (+) node!!!\n");
                    ovp_coverage_push_node(pos_g_idx, pos_g_num, pos_n_idx, pos_n_num, pos_reach_end);
                    // Cov of the positive node
                    segs_overlap_chain('+', pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend,bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, end_chaining_cov, bundle->bpcov);
                    target_cov_ratio = end_chaining_cov/expected_cov_pos;
                    fprintf(stderr, "$$$$ Ratio (+): %f\n", target_cov_ratio);
                } else {
                    // Positive reaches the end. Push negative all the way to the end
                    if (pos_reach_end && !neg_reach_end) {
                        fprintf(stderr, "\t>>> Pushing (-) node!!!\n");
                        ovp_coverage_push_node(neg_g_idx, neg_g_num, neg_n_idx, neg_n_num, neg_reach_end);
                        // Cov of the negative node
                        segs_overlap_chain('-', pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend,bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, end_chaining_cov, bundle->bpcov);
                        target_cov_ratio = end_chaining_cov/expected_cov_neg;
                        fprintf(stderr, "$$$$ Ratio (-): %f\n", target_cov_ratio);
                    } else if (pos_reach_end && neg_reach_end) {
                        // Shouldn't enter this chunk.
                    }
                }
            } else if (pos_n_end > neg_n_end) {
                // Push negative node
                if (!neg_reach_end) {
                    fprintf(stderr, "\t>>> Pushing (-) node!!!\n");
            // fprintf(stderr, "Before ?>> neg_g_idx: %d; neg_g_num: %d; neg_n_idx: %d; neg_n_num: %d; neg_reach_end: %d\n", neg_g_idx, neg_g_num, neg_n_idx, neg_n_num, neg_reach_end);

            // fprintf(stderr, "Before ?>> pos_g_idx: %d; pos_g_num: %d; pos_n_idx: %d; pos_n_num: %d; pos_reach_end: %d\n", pos_g_idx, pos_g_num, pos_n_idx, pos_n_num, pos_reach_end);
                    ovp_coverage_push_node(neg_g_idx, neg_g_num, neg_n_idx, neg_n_num, neg_reach_end);

            // fprintf(stderr, "After ?>> neg_g_idx: %d; neg_g_num: %d; neg_n_idx: %d; neg_n_num: %d; neg_reach_end: %d\n", neg_g_idx, neg_g_num, neg_n_idx, neg_n_num, neg_reach_end);

            // fprintf(stderr, "After ?>> pos_g_idx: %d; pos_g_num: %d; pos_n_idx: %d; pos_n_num: %d; pos_reach_end: %d\n", pos_g_idx, pos_g_num, pos_n_idx, pos_n_num, pos_reach_end);

                    // Cov of the negative node
                    segs_overlap_chain('-', pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend,bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, end_chaining_cov, bundle->bpcov);
                    target_cov_ratio = end_chaining_cov/expected_cov_neg;
                    fprintf(stderr, "$$$$ Ratio (-): %f\n", target_cov_ratio);
                } else {
                    // Negative reaches the end. Push positive all the way to the end
                    fprintf(stderr, "Negative reaches the end. Push positive all the way to the end!!\n");
                    if (neg_reach_end && !pos_reach_end) {
                        fprintf(stderr, "\ts>>> Pushing (+) node!!!\n");
                        ovp_coverage_push_node(pos_g_idx, pos_g_num, pos_n_idx, pos_n_num, pos_reach_end);
                        // Cov of the positive node
                        segs_overlap_chain('+', pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend,bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, end_chaining_cov, bundle->bpcov);
                        target_cov_ratio = end_chaining_cov/expected_cov_pos;
                        fprintf(stderr, "$$$$ Ratio (+): %f\n", target_cov_ratio);
                    } else if (pos_reach_end && neg_reach_end) {
                        // Shouldn't enter this chunk.
                    }
                }
            } else if (pos_n_end == neg_n_end) {
                // Only push one node => either positive/negative
                // Push positive first!!!
                if (!pos_reach_end) {
                    fprintf(stderr, ">>> Pushing (+) node!!!\n");
                    ovp_coverage_push_node(pos_g_idx, pos_g_num, pos_n_idx, pos_n_num, pos_reach_end);
                    // Cov of the positive node
                    segs_overlap_chain('+', pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend,bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, end_chaining_cov, bundle->bpcov);
                    target_cov_ratio = end_chaining_cov/expected_cov_pos;
                    fprintf(stderr, "$$$$ Ratio (+): %f\n", target_cov_ratio);
                } else {
                    if (!neg_reach_end) {
                        fprintf(stderr, ">>> Pushing (-) node!!!\n");
                        ovp_coverage_push_node(neg_g_idx, neg_g_num, neg_n_idx, neg_n_num, neg_reach_end);
                        // Cov of the negative node
                        segs_overlap_chain('-', pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend,bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, last_ovp_end, chaining_hold_strand, chaining_hold_cov, end_chaining_cov, bundle->bpcov);
                        target_cov_ratio = end_chaining_cov/expected_cov_neg;
                        fprintf(stderr, "$$$$ Ratio (-): %f\n", target_cov_ratio);
                    }
                }
            }

            if (target_cov_ratio > 2 || target_cov_ratio < 0.5) {
                fprintf(stderr, "$$$$ Warning!!!!!! The target_cov_ratio is out of range!!!!\n");
            }
            fprintf(stderr, "\n\n");
            ratio_file << target_cov_ratio << "\n";
        }
        






        // while (pos_g_idx < pos_g_num || neg_g_idx < neg_g_num) {
        //     // if ((pos_g_idx >= graphs_vec[1]->Count()) || (neg_g_idx >= graphs_vec[0]->Count())) {

        //     //     break;
        //     //     if ((pos_g_idx >= graphs_vec[1]->Count()) && (neg_g_idx >= graphs_vec[0]->Count())) {
        //     //         break;
        //     //     } else {
        //     //         if (pos_g_idx >= graphs_vec[1]->Count()) {
        //     //             pos_g_idx = graphs_vec[1]->Count() - 1;
        //     //             // neg_g_idx++;
        //     //             // neg_n_idx = 1;
        //     //         }
        //     //         if (neg_g_idx >= graphs_vec[0]->Count()) {
        //     //             neg_g_idx = graphs_vec[0]->Count()-1;
        //     //             // pos_g_idx++;
        //     //             // pos_n_idx = 1;
        //     //         }
        //     //     }
        //     // }

        //     if ((pos_g_idx == pos_g_num-1) || (neg_g_idx == neg_g_num-1)) {
        //         break;
        //     }
        //     fprintf(stderr, "pos_g_idx: %d;  neg_g_idx: %d\n", pos_g_idx, neg_g_idx);

        //     int pos_g_start = graphs_vec[1]->Get(pos_g_idx)->refstart;
        //     int pos_g_end = graphs_vec[1]->Get(pos_g_idx)->refend;

        //     int neg_g_start = graphs_vec[0]->Get(neg_g_idx)->refstart;
        //     int neg_g_end = graphs_vec[0]->Get(neg_g_idx)->refend;


        //     bool g_overlap = segs_overlap(pos_g_start, pos_g_end, neg_g_start, neg_g_end);
        //     if (g_overlap) {
        //         // Process the overlapping graphs

        //         fprintf(stderr, "graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1: %d;  graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1: %d\n", graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count(), graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count());

        //         // Use -1 for trick
        //         int pos_n_num = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Count()-1;
        //         int neg_n_num = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Count()-1;
        //         while (pos_n_idx < pos_n_num || neg_n_idx < neg_n_num) {

        //             if (pos_n_idx == pos_n_num-1 || neg_n_idx == neg_n_num-1) {
                        
        //                 // Positive graph reach the last node.
        //                 if (pos_n_idx < pos_n_num-1) {
        //                     // It means, 
        //                     // 1. pos_n_idx < pos_n_num-1
        //                     // 2. pos_g_idx ? pos_g_num-1
        //                     // 3. neg_n_idx == neg_n_num-1
        //                     // 4. neg_g_idx ? neg_g_num-1
        //                     pos_n_idx++;
        //                 } else if (pos_n_idx == pos_n_num-1) {

        //                     // This is the last positive graph in the graphs_vec
        //                     if (pos_g_idx == pos_g_num-1) {
        //                         if (neg_n_idx == neg_n_num-1) {
        //                             if (neg_g_idx == neg_g_num-1) {
        //                                 // It means, 
        //                                 // 1. pos_n_idx == pos_n_num-1
        //                                 // 2. pos_g_idx == pos_g_num-1
        //                                 // 3. neg_n_idx == neg_n_num-1
        //                                 // 4. neg_g_idx == neg_g_num-1
        //                                 break;
        //                             } else {
        //                                 // It means, 
        //                                 // 1. pos_n_idx == pos_n_num-1
        //                                 // 2. pos_g_idx == pos_g_num-1
        //                                 // 3. neg_n_idx == neg_n_num-1
        //                                 // 4. neg_g_idx < neg_g_num-1
        //                                 neg_g_idx ++;
        //                                 neg_n_idx = 1;
        //                                 break;
        //                             }
        //                         } else if (neg_n_idx < neg_n_num-1) {
        //                                 // It means, 
        //                                 // 1. pos_n_idx == pos_n_num-1
        //                                 // 2. pos_g_idx == pos_g_num-1
        //                                 // 3. neg_n_idx < neg_n_num-1
        //                                 // 4. neg_g_idx < neg_g_num-1
        //                             neg_n_idx++;
        //                         } 
        //                     } else {
        //                         // It means, 
        //                         // 1. pos_n_idx == pos_n_num-1
        //                         // 2. pos_g_idx < pos_g_num-1
        //                         pos_g_idx++;
        //                         pos_n_idx = 1;
        //                         break;
        //                     }
        //                 }

        //                 // // Negative graph reach the last node.
        //                 // if (neg_n_idx == neg_n_num-1) {
        //                 //     // This is the last negative graph in the graphs_vec
        //                 //     if (neg_g_idx == neg_g_num-1) {
        //                 //         if (pos_n_idx == pos_n_num-1) {
        //                 //             if (pos_g_idx == pos_g_num-1) {
        //                 //                 // It means, 
        //                 //                 // 1. neg_n_idx == neg_n_num-1
        //                 //                 // 2. neg_g_idx == neg_g_num-1
        //                 //                 // 3. pos_n_idx == pos_n_num-1
        //                 //                 // 4. pos_g_idx == pos_g_num-1
        //                 //                 break;
        //                 //             } else {
        //                 //                 // It means, 
        //                 //                 // 1. neg_n_idx == neg_n_num-1
        //                 //                 // 2. neg_g_idx == neg_g_num-1
        //                 //                 // 3. pos_n_idx == pos_n_num-1
        //                 //                 // 4. pos_g_idx < pos_g_num-1
        //                 //                 pos_g_idx ++;
        //                 //                 pos_n_idx = 1;
        //                 //             }
        //                 //         } else if (pos_n_idx < pos_n_num-1) {
        //                 //                 // It means, 
        //                 //                 // 1. neg_n_idx == neg_n_num-1
        //                 //                 // 2. neg_g_idx == neg_g_num-1
        //                 //                 // 3. pos_n_idx < pos_n_num-1
        //                 //                 // 4. pos_g_idx < pos_g_num-1
        //                 //             pos_n_idx++;
        //                 //         } 
        //                 //     } else {
        //                 //         // It means, 
        //                 //         // 1. neg_n_idx == neg_n_num-1
        //                 //         // 2. neg_g_idx < neg_g_num-1
        //                 //         neg_g_idx++;
        //                 //         neg_n_idx = 1;
        //                 //         break;
        //                 //     }
        //                 // }
        //             }

        //             uint pos_n_start = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->start;
        //             uint pos_n_end = graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->end;
        //             uint neg_n_start = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->start;
        //             uint neg_n_end = graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->end;
        //             bool n_overlap = segs_overlap(pos_n_start, pos_n_end, neg_n_start, neg_n_end);
        //             if (n_overlap) {
        //                 // Process the overlapping nodes


        //                 // uint overlap_n_start = pos_n_start<neg_n_start ? pos_n_start : neg_n_start;
        //                 // uint overlap_n_end = pos_n_end>neg_n_end ? pos_n_end : neg_n_end;
        //                 // This is the the unstranded reads coverage
        //                 // node_coverage_neg = get_cov(0, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
        //                 node_coverage_uns = get_cov(1, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
        //                 // node_coverage_pos = get_cov(2, pos_n_start-refstart, pos_n_end-refstart, bundle->bpcov);
        //                 float node_coverage_pos = 0.0;
        //                 float node_coverage_neg = 0.0;

        //                 calculate_ovp_coverage(pos_n_start, pos_n_end, neg_n_start, neg_n_end, refstart, refend, bundle_coverage_ratio_pos, bundle_coverage_ratio_neg, node_coverage_pos, node_coverage_neg, bundle->bpcov);

        //                 // node_coverage = get_cov(1, overlap_n_start-refstart, overlap_n_end-refstart, bundle->bpcov);

        //                 fprintf(stderr, "pos_g_idx: %d; pos_n_idx: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\nneg_g_idx: %d; neg_n_idx: %d; nodeid: %d (%u - %u) => node_coverage (%f, %f / %f, %f)\nunstranded node_coverage (%f)\n", pos_g_idx, pos_n_idx, graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->nodeid, pos_n_start, pos_n_end, node_coverage_pos, node_coverage_pos/(pos_n_end-pos_n_start), graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->cov_s->Get(0), graphs_vec[1]->Get(pos_g_idx)->no2gnode_unispg[1][0].Get(pos_n_idx)->cov_s->Get(0)/(pos_n_end-pos_n_start), neg_g_idx, neg_n_idx, graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->nodeid, neg_n_start, neg_n_end, node_coverage_neg, node_coverage_neg/(neg_n_end-neg_n_start), graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->cov_s->Get(0), graphs_vec[0]->Get(neg_g_idx)->no2gnode_unispg[0][0].Get(neg_n_idx)->cov_s->Get(0)/(neg_n_end-neg_n_start), node_coverage_uns);

        //                 if (pos_n_end < neg_n_end) {
        //                     // positive ends first
        //                     pos_n_idx++;
        //                 } else if (neg_n_end < pos_n_end) {
        //                     // negative ends first
        //                     neg_n_idx++;
        //                 } else if (neg_n_end == pos_n_end) {
        //                     pos_n_idx++;
        //                     neg_n_idx++;
        //                 }


        //             } else {
        //                 if (pos_n_end < neg_n_start) {
        //                     // Positive is in the front & not overlap
        //                     pos_n_idx++;               
        //                     if (pos_n_idx >= pos_g_num-1) {
        //                         // This is the last positive graph in the graphs_vec
        //                         if (neg_g_idx < neg_g_num-1) {
        //                             neg_g_idx++;
        //                         }
        //                     } else {
        //                         pos_n_idx++;               
        //                     }




        //                 } else if (neg_n_end < pos_n_start) {
        //                     // negative is in the front & not overlap
        //                     neg_n_idx++;
        //                 }
        //             }
        //         }







        //         if (pos_g_end < neg_g_end) {
        //             // positive ends first
        //             pos_g_idx++;
        //         } else if (neg_g_end < pos_g_end) {
        //             // negative ends first
        //             neg_g_idx++;
        //         } else if (neg_g_end == pos_g_end) {
        //             pos_g_idx++;
        //             neg_g_idx++;
        //         }






                
        //     } else {
        //         if (pos_g_end < neg_g_start) {
        //             // Positive is in the front & not overlap
        //             if (pos_g_idx >= pos_g_num-1) {
        //                 // This is the last positive graph in the graphs_vec
        //                 if (neg_g_idx < neg_g_num-1) {
        //                     neg_g_idx++;
        //                 }
        //             } else {
        //                 pos_g_idx++;               
        //             }
        //         } else if (neg_g_end < pos_g_start) {
        //             // negative is in the front & not overlap
        //             if (neg_g_idx >= neg_g_num-1) {
        //                 // This is the last negative graph in the graphs_vec
        //                 if (pos_g_idx < pos_g_num-1) {
        //                     pos_g_idx++;
        //                 }
        //             } else {
        //                 neg_g_idx++;
        //             }
        //         }
        //     }
        // }
    
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
