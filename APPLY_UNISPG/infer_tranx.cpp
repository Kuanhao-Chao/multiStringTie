#include "infer_tranx.h"

void infer_transcripts_unispg(BundleData* bundle, UniSpliceGraph* unispg) {

    count_good_junctions(bundle);
	
    // geneno = build_graphs_unispg(bundle, unispg);


    fprintf(stderr, "Inside `infer_transcripts_unispg`\n");
    fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());


    for (int i=0; i<unispg->get_no2gnode()->Count(); i++) {
        fprintf(stderr, "unispg[%d] nodeid: %d\n", i, unispg->get_no2gnode()->Get(i)->nodeid);
    }


    for (int n=0;n<bundle->readlist.Count();n++) {
        fprintf(stderr, "bundle->readlist[%d]: start:%u - end:%u \n ", n , bundle->readlist[n]->start, bundle->readlist[n]->end);
        fprintf(stderr, "n: %d; bundle->bpcov[0].Count(): %d  , bundle->bpcov[1].Count(): %d  , bundle->bpcov[2].Count(): %d  \n ", n , bundle->bpcov[0].Count(), bundle->bpcov[1].Count(), bundle->bpcov[2].Count());

        float single_count=bundle->readlist[n]->read_count;

        // for(int j=0; j<bundle->readlist[n]->pair_idx.Count();j++) {
        //     int np=bundle->readlist[n]->pair_idx[j];

        //     fprintf(stderr, "n: %d  np: %d ", n, np);

        //     // if(np>-1) {
        //     //     single_count-=bundle->readlist[n]->pair_count[j];
        //     //     if(n<np) {
        //     //         get_fragment_pattern(bundle->readlist,n,np,bundle->readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
        //     //     }
        //     // }
        // }
        // if(single_count>epsilon) {
        //     get_fragment_pattern(bundle->readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
        // }

    }
}
