#include "infer_tranx_A.h"

void infer_transcripts_APPLY_UNISPG(BundleData* bundle, GPVec<UnispgGp_APPLY>** graphs_vec) {

    int refstart = bundle->start;
    int refend = bundle->end;

    count_good_junctions_APPLY_UNISPG(bundle);
	
    // geneno = build_graphs_unispg(bundle, unispg);
    fprintf(stderr, "Inside `infer_transcripts_APPLY_UNISPG`\n");
    fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
    fprintf(stderr, "bundle->bpcov[1].Count(): %d\n", bundle->bpcov[1].Count());
    fprintf(stderr, "bundle->start - end: %d - %d\n", bundle->start, bundle->end);


    /*****************************
     * Iterate through graphs.
     *****************************/
    for(int s=0;s<2;s+=1) {
        fprintf(stderr, "\t@@ graphs_vec[%d]: %d\n", s, graphs_vec[s]->Count());
        for (int i=0; i<graphs_vec[s]->Count(); i++) {
            fprintf(stderr, "\t\t graphs_vec bound: %u - %u \n", graphs_vec[s]->Get(i)->get_refstart(), graphs_vec[s]->Get(i)->get_refend());


        }
    }

    /*****************************
     **    'get_fragment_pattern_APPLY_UNISPG' function
     **        because of this going throu
     *****************************/
    int global_gidx[2] = {0};
    for (int n=0; n < bundle->readlist.Count(); n++) {
    	float single_count=bundle->readlist[n]->read_count;

        uint read_start = bundle->readlist[n]->start;
        uint read_end = bundle->readlist[n]->end;

        // Only process the reads that are completely in the graph boundary.
        if (read_start >= refstart && read_end <= refend) {
            for (int i=0; i<bundle->readlist[n]->pair_idx.Count(); i++) {
                // fprintf(stderr, ">> single_count: %d\n", i);
                int np = bundle->readlist[n]->pair_idx[i];
                if (np > -1) {
                    single_count -= bundle->readlist[n]->pair_idx[i];
                    fprintf(stderr, ">> n: %d;   np: %d\n", n, np);
                    if (n < np) {
                        // fprintf(stderr, ">> n < np: %d\n", n < np);
                        get_fragment_pattern_APPLY_UNISPG(bundle, bundle->readlist, n, np, bundle->readlist[n]->pair_count[i], graphs_vec, global_gidx);
                            // readlist,n,np,readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
                    }
                }
            }
            if (single_count > epsilon) {
                get_fragment_pattern_APPLY_UNISPG(bundle, bundle->readlist, n, -1, single_count, graphs_vec, global_gidx);
                    // readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
            }
        } else {
            skip_counter += 1;
            // bundle->readlist[n]->juncs;
            fprintf(stderr, ">> nh: %d\n", bundle->readlist[n]->nh);
            // bundle->readlist[n]->pair_count;
            // bundle->readlist[n]->pair_idx;
            fprintf(stderr, ">> read_count: %f\n", bundle->readlist[n]->read_count);
            // bundle->readlist[n]->segs;
            fprintf(stderr, ">> strand: %d\n", bundle->readlist[n]->strand);

            fprintf(stderr, "Skipped!!!\n");
        }
    }



    
    

    // /*****************************
    //  ** Output the Stringtie vs redistributed coverage 
    //  *****************************/   
    if (graphs_vec[0]->Count() == 0) {
        fprintf(stderr, "####### Positive strand!!!\n");
        for (int g=0; g<graphs_vec[1]->Count(); g++) {
            for (int n=1; n<graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Count()-1; n++) {
                float target_cov_ratio = 0.0;
                fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[1]->Get(g)->get_refstart(), graphs_vec[1]->Get(g)->get_refend(), g, n, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->nodeid, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->start, graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->end);

                

                float expected_cov_pos = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_s->Get(0);
                int expected_len_pos = int(graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->len());
                float new_cov_pos = graphs_vec[1]->Get(g)->no2gnode_unispg[1][0].Get(n)->cov_unispg_s[0];

                fprintf(stderr, ">> expected_cov_pos: %f;  new_cov_pos: %f\n", expected_cov_pos/expected_len_pos, new_cov_pos/expected_len_pos);

                cov_file_pos << expected_cov_pos << "\t" << new_cov_pos << "\n";
                cov_file_pos_norm << (expected_cov_pos/expected_len_pos) << "\t" << (new_cov_pos/expected_len_pos) << "\n";
            }
        }
    } else if (graphs_vec[1]->Count() == 0) {
        fprintf(stderr, "####### Negative strand!!!\n");
        for (int g=0; g<graphs_vec[0]->Count(); g++) {
            for (int n=1; n<graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Count()-1; n++) {
                float target_cov_ratio = 0.0;
                fprintf(stderr, "g: %d; n: %d; ref (%d - %d); graphs_vec[%d][%d] nodeid: %d (%u - %u)\n", g, n, graphs_vec[0]->Get(g)->get_refstart(), graphs_vec[0]->Get(g)->get_refend(), g, n, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->nodeid, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->start, graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->end);
                
                float expected_cov_neg = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_s->Get(0);
                int expected_len_neg = int(graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->len());
                float new_cov_neg = graphs_vec[0]->Get(g)->no2gnode_unispg[0][0].Get(n)->cov_unispg_s[0];

                fprintf(stderr, ">> expected_cov_neg: %f;  new_cov_neg: %f.", expected_cov_neg/expected_len_neg, new_cov_neg/expected_len_neg);
                

                cov_file_neg << expected_cov_neg << "\t" << new_cov_neg << "\n";
                cov_file_neg_norm << (expected_cov_neg/expected_len_neg) << "\t" << (new_cov_neg/expected_len_neg) << "\n";
            }
        }
    }
}
