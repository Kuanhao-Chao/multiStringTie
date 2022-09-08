#include "unispg.h"
#include "GBitVec.h"
#include <float.h>
#include <limits.h>

#include <iostream>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

extern UnispgGp* unispg_gp;
extern GVec<int> current_gidx;
extern FILE* uinigraph_out;

void printBitVecTest(GBitVec& bv) {
   for (uint i=0;i<bv.size();i++) {
       fprintf(stderr, "%c", bv.test(i)?'1':'0');
   }
}

void UnispgGp::ProcessSample(GStr sample_name) {
    samples.Add(sample_name);
    // for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
    //     int s=sno/2; // adjusted strand due to ignoring neutral strand
    //     current_gidx[s] = 0;
    //     new_gidx[s] = 0;
    //     last_nidx[s] = 1; // node id
    //     prev_bdy[s] = 0;
    //     lclg_bundle_num[s] = 1;
    //     has_unispg_tail[s] = false;

    //     // source_gp[s] = new CGraphnodeUnispg();
    //     // sink_gp[s] = new CGraphnodeUnispg();
    // }
}

// void UnispgGp::PrintGraphGp() {
//     fprintf(stderr, "*********************************\n");
//     fprintf(stderr, "*********** PrintGraphGp ********\n");
//     fprintf(stderr, "*********************************\n");

//     { // sDEBUG ONLY
//         printTime(stderr);
//         // for(int s=0;s<2;s++) {
//         //     fprintf(stderr, "\n\tThere are %d stranded[%d] graphs\n", gpSize[s],int(2*s));
//         //     // if () {

//         //     // }
//         //     for(int b=0;b<gpSize[s];b++) {
//         //         fprintf(stderr, ">>>>>>> 1-2 gpSize[%d]: %d\n", s, gpSize[s]);
//         //         fprintf(stderr, ">>>>>>> 1-2 graphnoGp[%d][%d]: %d\n", s, b, graphnoGp[s][b]);
//         //         if(graphnoGp[s][b]) {
//         //             GStr pat;
//         //             fprintf(stderr,"\t\tGraph[%d][%d] with %d nodes and %d edges :",int(2*s),b,graphnoGp[s][b],edgenoGp[s][b]);
//         //             for(int nd=1;nd<graphnoGp[s][b]-1;nd++) {
//         //                 fprintf(stderr, "&&& nd: %d\n", nd);
//         //                 fprintf(stderr, "&&& no2gnodeGp[s][b][nd]->start: %d\n", no2gnodeGp[s][b][nd]->start);
//         //                 fprintf(stderr, "&&& no2gnodeGp[s][b][nd]->end: %d\n", no2gnodeGp[s][b][nd]->end);
//         //                 fprintf(stderr," %d(%d-%d)",nd,no2gnodeGp[s][b][nd]->start,no2gnodeGp[s][b][nd]->end);
//         //             }
//         //         } 
//         //         fprintf(stderr,"\n");
//         //     }
//         // }

//     //     for(int s=0;s<2;s++) {
//     // 	fprintf(stderr, "There are %d stranded[%d] graphs\n",bno[s],int(2*s));
//     // 	// fprintf(stderr, "1 bundle[sno].Count(): %d\n", bno[s]);
//     // 	for(int b=0;b<bno[s];b++) {
//     // 		if(graphno[s][b]) {
//     // 			GStr pat;
//     // 			fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges with lastgpos=%d:",int(2*s),b,graphno[s][b],edgeno[s][b],lastgpos[s][b]);
//     // 			for(int nd=1;nd<graphno[s][b]-1;nd++)
//     // 				fprintf(stderr," %d(%d-%d)",nd,no2gnode[b][nd]->start,no2gnode[b][nd]->end);
//     // 			fprintf(stderr,"\n");
//     // 			print_pattern(tr2no[s][b],pat,graphno[s][b]);
//     // 		}
//     // 	}
//     // }
//     }
// }

GPVec<CGraphnodeUnispg>** UnispgGp::get_no2gnodeGp () {
    return no2gnode_unispg;
}