#include "unispg.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif



void unispgGp_clean() {
    for(int i=0;i<2;i++) {
    // gpSize[i] = 0;
    // graphnoGp[i].Clear();
    // edgenoGp[i].Clear();
    delete [] no2gnodeGp_unispg[i];
    // no2gnodeGp[i] = new GPVec<CGraphnode>;
    // no2gnode[i] = NULL;
    };
}

// void SetUnispgCapacity(int s) {
//     fprintf(stderr, "**************************************\n");
//     fprintf(stderr, "*********** SetUnispgCapacity ********\n");
//     fprintf(stderr, "**************************************\n");
//     // fprintf(stderr, "s: %d; capacity: %d\n", s, capacity);
//     // gpSize[s] = 2;
//     // graphnoGp[s].setCapacity(capacity);
//     // edgenoGp[s].setCapacity(capacity);
//     // graphnoGp[s].cAdd(0);
//     // edgenoGp[s].cAdd(0);
//     no2gnodeGp[s] = new GPVec<CGraphnode>[20000];
// }

// void SetRefStartEnd(int refstart_i, int refend_i) {
//     refstart = refstart_i;
//     refend = refend_i;
// }

// void UpdateRefStartEnd(int refstart_i, int refend_i) {
//     refstart = refstart_i;
//     refend = refend_i;
// }

void AddGraph(int fidx, int s, int track_idx, GPVec<CGraphnode>* no2gnode) {
    int cgidx = current_gidx[s];
    fprintf(stderr, "*****************************\n");
    fprintf(stderr, "*********** AddGraph ********\n");
    fprintf(stderr, "*****************************\n");
        // GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
        // GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
        // GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i


    // fprintf(stderr, "* uni_splice_graph.refstart: %d \n", uni_splice_graph -> get_refstart());
    // fprintf(stderr, "* uni_splice_graph.refend: %d \n", uni_splice_graph -> get_refend());
    // fprintf(stderr, "* uni_splice_graph.s: %d \n", uni_splice_graph -> get_s());
    // fprintf(stderr, "* uni_splice_graph.g_idx: %d \n", uni_splice_graph -> get_g_idx());
    // fprintf(stderr, "* uni_splice_graph.get_graphno: %d \n", uni_splice_graph -> get_graphno());
    // fprintf(stderr, "* uni_splice_graph.get_edgeno: %d \n", uni_splice_graph -> get_edgeno());

    // int refstart;
    // int refend;
    // int s;
    // int g_idx;

    // int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    // int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    // GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i

    for (int i=0; i<no2gnode->Count(); i++) {
        CGraphnode* node = new CGraphnode(no2gnode->Get(i));
        // fprintf(stderr, "New node address %p\n", node);
        // fprintf(stderr, "New node sizeof  %d\n", sizeof(node));
        // fprintf(stderr, "Old node address %p\n", no2gnode->Get(i));
        // fprintf(stderr, "Old node sizeof %p\n", sizeof(no2gnode->Get(i)));

        node = no2gnode->Get(i);
        // fprintf(stderr, "Before no2gnodeGp_unispg[s][0] graph address %p\n", &(no2gnodeGp_unispg[s][0]));
        // fprintf(stderr, "Before no2gnodeGp_unispg[s][0] sizeof  %d\n", sizeof(no2gnodeGp_unispg[s][0]));
        // fprintf(stderr, "Before no2gnodeGp_unispg[s][0].Count()  %d\n", no2gnodeGp_unispg[s][0].Count());
        no2gnodeGp_unispg[s][cgidx].Add(node);
        // fprintf(stderr, "After no2gnodeGp_unispg[s][0] graph address %p\n", &no2gnodeGp_unispg[s][0]);
        // fprintf(stderr, "After no2gnodeGp_unispg[s][0] sizeof  %d\n", sizeof(no2gnodeGp_unispg[s][0]));
        // fprintf(stderr, "After no2gnodeGp_unispg[s][0].Count()  %d\n", no2gnodeGp_unispg[s][0].Count());

        // fprintf(stderr, "Inside!! cgidx: %d \n", cgidx);
        // fprintf(stderr, "no2gnode->Get(i)->nodeid: %d \n", node->nodeid);
        // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx].Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx].Last()->nodeid);

        // fprintf(stderr, "no2gnodeGp_unispg[s][0] graph address: %p \n", &(no2gnodeGp_unispg[s][0]));
        
        // fprintf(stderr, "no2gnodeGp_unispg[s][0] sizeof: %d \n", sizeof(no2gnodeGp_unispg[s][0]));

        // fprintf(stderr, "no2gnodeGp_unispg[s][0].Count(): %d \n", no2gnodeGp_unispg[s][0].Count());
        // fprintf(stderr, "no2gnodeGp_unispg[s][0]->Last()->nodeid: %d \n", no2gnodeGp_unispg[s][0].Last()->nodeid);


        // fprintf(stderr, "Inside!! cgidx: %d \n", cgidx);
        // fprintf(stderr, "no2gnode->Get(i)->nodeid: %d \n", node->nodeid);

        // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx].Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx].Last()->nodeid);
        // fprintf(stderr, "no2gnodeGp_unispg[s][0] graph address: %p \n", &(no2gnodeGp_unispg[s][0]));
        
        // fprintf(stderr, "no2gnodeGp_unispg[s][0] sizeof: %d \n", sizeof(no2gnodeGp_unispg[s][0]));

        // fprintf(stderr, "no2gnodeGp_unispg[s][0].Count(): %d \n", no2gnodeGp_unispg[s][0].Count());
        // fprintf(stderr, "no2gnodeGp_unispg[s][0]->Last()->nodeid: %d \n", no2gnodeGp_unispg[s][0].Last()->nodeid);

        // if (cgidx >= 1) {
        //     fprintf(stderr, "no2gnodeGp_unispg[s][cgidx-1]->Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx-1].Last()->nodeid);
        //     fprintf(stderr, "no2gnodeGp_unispg[s][cgidx]->Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx].Last()->nodeid);
        //     // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx+1]->Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx+1].Last()->nodeid);
        // }
    }

    if (fidx == 0) {
        // This is the first unispg. Simply add copy it into UnispgGp
        // graphnoGp[s].Add(graphno);
        // edgenoGp[s].Add(edgeno);
        // gpSize[s] = no2gnode->Count();
        // no2gnodeGp_unispg[s][b].Add(no2gnode->Get(0));

        // for (int i=0; i<no2gnode->Count(); i++) {
        //     no2gnodeGp_unispg[s][cgidx].Add(no2gnode->Get(i));
        //     // cgidx = cgidx+1;
        //     fprintf(stderr, "Inside!! cgidx: %d \n", cgidx);
        // }   
        
        // fprintf(stderr, "$$$      no2gnodeGp_unispg[s]: %d\n", sizeof(no2gnodeGp_unispg[s]));

        // for (int i = 0; i < no2gnode->Count(); i++) {
        //     // no2gnode->Get(i);
        //     fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
        // }

        // for (int i = 0; i < no2gnodeGp_unispg[s][cgidx].Count(); i++) {
        //     fprintf(stderr, "&&&& This is the global graphnode: %d\n", no2gnodeGp_unispg[s][cgidx][i]->nodeid);
        // }

    } else {
        fprintf(stderr, "************************************\n");
        fprintf(stderr, "*********** Add more graph! ********\n");
        fprintf(stderr, "************************************\n");
        // Need to check overlapping & new graph creation

        // for(int sno=0; sno<3; sno+=2) { // skip neutral bundles -> those shouldn't have junctions
        //     int s=sno/2; // adjusted strand due to ignoring neutral strand
        //     // fprintf(stderr, "2 bundle[sno].Count(): %d\n", bundle[sno].Count());
        //     for(int b=0;b<bundle[sno].Count();b++) {
        //     }
        // }



        // for (int i = 0; i < no2gnode->Count(); i++) {
        //     // no2gnode->Get(i);
        //     fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
        // }
        // for (int i = 0; i < no2gnodeGp_unispg[s][cgidx].Count(); i++) {
        //     fprintf(stderr, "&&&& This is the global graphnode: %d\n", no2gnodeGp_unispg[s][cgidx][i]->nodeid);
        // }
        // track_idx++;


        int lclg_start = no2gnode->Get(1)->start;
        int lclg_end = no2gnode->Get(no2gnode->Count()-2)->end;

        int unispg_start = no2gnodeGp_unispg[s][cgidx][1]->start;
        int unispg_end = no2gnodeGp_unispg[s][cgidx][ no2gnodeGp_unispg[s][cgidx].Count()-2 ]->end;

        fprintf(stderr, "$$$      cgidx: %d\n", cgidx);
        fprintf(stderr, "$$$      no2gnodeGp_unispg[s]: %d\n", sizeof(no2gnodeGp_unispg[s]));
        fprintf(stderr, "$$$  track_idx: %d\n", track_idx);
        fprintf(stderr, "$$$ lclg_start: %d,  lclg_end: %d,  unispg_start: %d,  unispg_end: %d\n", lclg_start, lclg_end, no2gnodeGp_unispg[s][cgidx][1]->start, no2gnodeGp_unispg[s][cgidx][ no2gnodeGp_unispg[s][cgidx].Count()-2 ]->end);

        track_idx++;


        if (unispg_end < lclg_start) {
            // ----------   |(s).................(e)|
            fprintf(stderr,"\n  &&& Bundle: ----------   |(s).................(e)| \n");
            track_idx++;
        } else if (unispg_start < lclg_start && unispg_end >= lclg_start && unispg_end <= lclg_end) {
            // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
            fprintf(stderr,"\n  &&& Bundle: ----------|(s).................(e)|   or   -----|(s)-----............(e)| \n");
        } else if (unispg_start < lclg_start && unispg_end > lclg_end) {
            // -----|(s)------------(e)|--
            fprintf(stderr,"\n &&& Bundle: -----|(s)------------(e)|-- \n");
        } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
            // |(s)----------.................(e)| 
            fprintf(stderr,"\n &&& Bundle: |(s)----------............(e)|\n");
        } else if (unispg_start > lclg_start && unispg_end < lclg_end) {
            // |(s)........----------........(e)|
            fprintf(stderr,"\n |(s)........----------........(e)| \n");
        } else if (unispg_start > lclg_start && unispg_end == lclg_end) {
            // |(s)............----------(e)|
            fprintf(stderr,"\n &&& Bundle: |(s)............----------(e)| \n");
        } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
            // |(s)----------(e)|
            fprintf(stderr,"\n &&& Bundle: |(s)----------(e)| \n");
        } else if (unispg_start <= lclg_end && unispg_end > lclg_end) {
            // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
            fprintf(stderr,"\n &&& Bundle: (s)...............------(e)|-----    or   |(s).................(e)|---------- \n");
        } else if (unispg_start > lclg_end) {
            // The node is outside the current bundle => This node belongs to the next bundlenode
            // |(s).................(e)|   ----------
            fprintf(stderr,"\n &&& Bundle: |(s).................(e)|   ---------- \n");
        } else {
            fprintf(stderr,"\n &&& Unknown area!!!! \n");
        }

        for (int i = 0; i < no2gnode->Count(); i++) {
            // no2gnode->Get(i);
            fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
        }

        for (int i = 0; i < no2gnodeGp_unispg[s][10000].Count(); i++) {
            fprintf(stderr, "&&&& This is the global graphnode: %d\n", no2gnodeGp_unispg[s][10000][i]->nodeid);
        }
        // This is the local graph.
        *no2gnode;
        // This is the global graph.
        no2gnodeGp_unispg[s][track_idx];
    }

    // int g_idx = 0;
    // for (int i=0; i<no2gnode->Count(); i++) {
    //     no2gnodeGp_unispg[s][g_idx].Add(no2gnode->Get(i));
    //     g_idx++;
    // }


    // gpSize[uni_splice_graph->get_s()] = uni_splice_graph->get_g_idx()+1;
    // fprintf(stderr, "&&& uni_splice_graph->get_g_idx()+1: %d\n", uni_splice_graph->get_g_idx()+1);
    // fprintf(stderr, "&&& gpSize[uni_splice_graph->get_s()]: %d\n", gpSize[uni_splice_graph->get_s()]);
    // We need to reset graph_idx when (1)the new strands     
}


void Clear() {
    // fprintf(stderr, "**** Start Clearing !!!! \n ");
    for(int i=0;i<2;i++) {
    // gpSize[i] = 0;
    // graphnoGp[i].Clear();
    // graphnoGp[i].setCapacity(8192);
    // edgenoGp[i].Clear();
    // edgenoGp[i].setCapacity(8192);
    delete [] no2gnodeGp_unispg[i];
    no2gnodeGp_unispg[i] = new GPVec<CGraphnode>[8192];
    // no2gnodeGp_unispg[i]->setCapacity(8192);
    };
}

// Need to be modifed!!!
void PrintGraphGp() {
    fprintf(stderr, "*********************************\n");
    fprintf(stderr, "*********** PrintGraphGp ********\n");
    fprintf(stderr, "*********************************\n");

    { // DEBUG ONLY
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
        //                 fprintf(stderr, "&&& no2gnodeGp_unispg[s][b][nd]->start: %d\n", no2gnodeGp_unispg[s][b][nd]->start);
        //                 fprintf(stderr, "&&& no2gnodeGp_unispg[s][b][nd]->end: %d\n", no2gnodeGp_unispg[s][b][nd]->end);
        //                 fprintf(stderr," %d(%d-%d)",nd,no2gnodeGp_unispg[s][b][nd]->start,no2gnodeGp_unispg[s][b][nd]->end);
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
    // 				fprintf(stderr," %d(%d-%d)",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end);
    // 			fprintf(stderr,"\n");
    // 			print_pattern(tr2no[s][b],pat,graphno[s][b]);
    // 		}
    // 	}
    // }
    }
}