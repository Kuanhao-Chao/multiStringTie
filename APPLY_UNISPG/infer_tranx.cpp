#include "infer_tranx.h"

void infer_transcripts_unispg(BundleData* bundle, UnispgGp* unispg) {

    count_good_junctions(bundle);
	
    // geneno = build_graphs_unispg(bundle, unispg);
    fprintf(stderr, "Inside `infer_transcripts_unispg`\n");
    fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
    fprintf(stderr, "bundle->bpcov[1].Count(): %d\n", bundle->bpcov[1].Count());


    for (int i=0; i<unispg->no2gnode_unispg[1][0].Count(); i++) {
        fprintf(stderr, "unispg[%d] nodeid: %d\n", i, i<unispg->no2gnode_unispg[1][0].Get(i)->nodeid);
    }

	/*****************************
	 ** Applying reads to universal graphs.
	 *****************************/
    GList<CReadAln>& readlist = bundle->readlist;
    // GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
    // Create here!!
    GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    // Create here!!
    CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
    // Create here!!
    GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    // Create here!!
    GVec<int> lastgpos[2];
    // Input variable!!
    GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    // Input variable!!
    GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    // GVec<int> trnumber[2]; // how many transfrags are on a strand s, in a graph g -> I can find out this from transfrag[s][g].Count()
    // int ngraph[2]={0,0};   // how many graphs are in each strand: negative (0), or positive(1) -> keep one for each bundle


    // Input variable!!
    GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
    

    int bno[2]={0,0};
    
    
    // for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
        // int s=sno/2; // adjusted strand due to ignoring neutral strand
    int s = 1;
    //char strnd='-';
    //if(s) strnd='+';

    // bundle2graph[s]=NULL;
    // if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    transfrag[s]=NULL;
    no2gnode[s]=NULL;
    tr2no[s]=NULL;
    gpos[s]=NULL;

    uint unispg_nstart = 0;
    uint unispg_nend = 0;

    int refstart = unispg ->get_refstart();
    int refend = unispg->get_refend();

    int readlist_idx = 0;
    float node_coverage = 0.0;
    fprintf(stderr, "unispg->graphno_unispg[s]: %d\n", unispg->graphno_unispg[s]);
    fprintf(stderr, "unispg->edgeno_unispg[s]: %d\n", unispg->edgeno_unispg[s]);
    for (int i=1; i<unispg->graphno_unispg[s]-1; i++) {

        fprintf(stderr, "%d; ref (%d - %d); unispg[%d] nodeid: %d (%u - %u)\n", i, refstart, refend, i, unispg->no2gnode_unispg[1][0].Get(i)->nodeid, unispg->no2gnode_unispg[1][0].Get(i)->start, unispg->no2gnode_unispg[1][0].Get(i)->end);




        unispg_nstart = unispg->no2gnode_unispg[1][0].Get(i)->start;
        unispg_nend = unispg->no2gnode_unispg[1][0].Get(i)->end;

        // readlist[readlist_idx]->start;
        // readlist[readlist_idx]->end;
        
        node_coverage = get_cov(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);

        // node_coverage = get_cov_sign(0, unispg_nstart-refstart, unispg_nend-refstart, bundle->bpcov);
        fprintf(stderr, "unispg[%d] nodeid: %d (%u - %u) => node_coverage (%f / %f)\n", i, unispg->no2gnode_unispg[1][0].Get(i)->nodeid, unispg->no2gnode_unispg[1][0].Get(i)->start, unispg->no2gnode_unispg[1][0].Get(i)->end, node_coverage, unispg->no2gnode_unispg[1][0].Get(i)->cov_s->Get(0));
    }

    // }
    // for (int i=0; i<3; i++) {
    //     for (int pos=0; pos<bundle->end-bundle->start+1; pos++) {
    //         fprintf(stderr, "  >>> bpcov[%d][%d]: %f\n", i, pos, bundle->bpcov[i][pos]);
    //     }
    // }

    // for (int n=0;n<bundle->readlist.Count();n++) {
    //     fprintf(stderr, "bundle->readlist[%d]: start:%u - end:%u \n ", n , bundle->readlist[n]->start, bundle->readlist[n]->end);
    //     fprintf(stderr, "n: %d; bundle->bpcov[0].Count(): %d  , bundle->bpcov[1].Count(): %d  , bundle->bpcov[2].Count(): %d  \n ", n , bundle->bpcov[0].Count(), bundle->bpcov[1].Count(), bundle->bpcov[2].Count());

    //     float single_count=bundle->readlist[n]->read_count;

    //     // for(int j=0; j<bundle->readlist[n]->pair_idx.Count();j++) {
    //     //     int np=bundle->readlist[n]->pair_idx[j];

    //     //     fprintf(stderr, "n: %d  np: %d ", n, np);

    //     //     // if(np>-1) {
    //     //     //     single_count-=bundle->readlist[n]->pair_count[j];
    //     //     //     if(n<np) {
    //     //     //         get_fragment_pattern(bundle->readlist,n,np,bundle->readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
    //     //     //     }
    //     //     // }
    //     // }
    //     // if(single_count>epsilon) {
    //     //     get_fragment_pattern(bundle->readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
    //     // }

    // }
}
