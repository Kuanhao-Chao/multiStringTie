#include "unispg.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

extern UnispgGp* unispg_gp;
extern GVec<int> current_gidx;

// void unispgGp_clean() {
//     for(int i=0;i<2;i++) {
//     // gpSize[i] = 0;
//     // graphnoGp[i].Clear();
//     // edgenoGp[i].Clear();
//     // delete [] no2gnodeGp_unispg[i];
//     delete [] no2gnode_unispg[i];

//     // no2gnodeGp[i] = new GPVec<CGraphnode>;
//     // no2gnode[i] = NULL;
//     };
// }

// // void SetUnispgCapacity(int s) {
// //     fprintf(stderr, "**************************************\n");
// //     fprintf(stderr, "*********** SetUnispgCapacity ********\n");
// //     fprintf(stderr, "**************************************\n");
// //     // fprintf(stderr, "s: %d; capacity: %d\n", s, capacity);
// //     // gpSize[s] = 2;
// //     // graphnoGp[s].setCapacity(capacity);
// //     // edgenoGp[s].setCapacity(capacity);
// //     // graphnoGp[s].cAdd(0);
// //     // edgenoGp[s].cAdd(0);
// //     no2gnodeGp[s] = new GPVec<CGraphnode>[20000];
// // }

// // void SetRefStartEnd(int refstart_i, int refend_i) {
// //     refstart = refstart_i;
// //     refend = refend_i;
// // }

// // void UpdateRefStartEnd(int refstart_i, int refend_i) {
// //     refstart = refstart_i;
// //     refend = refend_i;
// // }

// void AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode) {
//         // GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
//         // GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
//         // GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i


//     // fprintf(stderr, "* uni_splice_graph.refstart: %d \n", uni_splice_graph -> get_refstart());
//     // fprintf(stderr, "* uni_splice_graph.refend: %d \n", uni_splice_graph -> get_refend());
//     // fprintf(stderr, "* uni_splice_graph.s: %d \n", uni_splice_graph -> get_s());
//     // fprintf(stderr, "* uni_splice_graph.g_idx: %d \n", uni_splice_graph -> get_g_idx());
//     // fprintf(stderr, "* uni_splice_graph.get_graphno: %d \n", uni_splice_graph -> get_graphno());
//     // fprintf(stderr, "* uni_splice_graph.get_edgeno: %d \n", uni_splice_graph -> get_edgeno());

//     // int refstart;
//     // int refend;
//     // int s;
//     // int g_idx;

//     // int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
//     // int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
//     // GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i

    

//     if (fidx == 0) {
//         fprintf(stderr, "*****************************\n");
//         fprintf(stderr, "*********** AddGraph ********\n");
//         fprintf(stderr, "*****************************\n");
//         int cgidx = current_gidx[s];
//         // This is the first unispg. Simply add copy it into UnispgGp
//         for (int i=0; i<no2gnode->Count(); i++) {
//             CGraphnode* node = new CGraphnode(no2gnode->Get(i));

//             GVec<bool> is_passed_s = NULL;
//             GVec<float> cov_s = NULL;
//             GVec<float> capacity_s = NULL;

//             CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(node->start, node->end, node->nodeid, is_passed_s, 1, cov_s, node->cov, capacity_s, node->capacity, node->rate);

//             // no2gnodeGp_unispg[s][cgidx].Add(node);
//             no2gnode_unispg[s][cgidx].Add(node_unispg);

//             // fprintf(stderr, "Inside!! cgidx: %d \n", cgidx);
//             // fprintf(stderr, "no2gnode->Get(i)->nodeid: %d \n", node->nodeid);
//             // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx].Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx].Last()->nodeid);
//             // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx].Last() start - end : %d - %d \n", int(no2gnodeGp_unispg[s][cgidx].Last()->start), int(no2gnodeGp_unispg[s][cgidx].Last()->end));
//             // fprintf(stderr, "no2gnode_unispg[s][cgidx].Last()->nodeid: %d \n", no2gnode_unispg[s][cgidx].Last()->nodeid);
//             // fprintf(stderr, "no2gnode_unispg[s][cgidx].Last() start - end : %d - %d \n", int(no2gnode_unispg[s][cgidx].Last()->start), int(no2gnode_unispg[s][cgidx].Last()->end));
//         }

//         // Linking parent and child
//         for (int i=0; i<no2gnode->Count(); i++) {
//             CGraphnode* node = new CGraphnode(no2gnode->Get(i));
//             fprintf(stderr, "CGraphnode parent: ");
// 			for(int p=0; p<no2gnode->Get(i)->parent.Count(); p++) {
//                 fprintf(stderr, "%d  ", no2gnode->Get(i)->parent[p]);
//                 fprintf(stderr, "%d  ", no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->parent[p] ]->nodeid);

//                 no2gnode_unispg[s][cgidx][i]->parent.Add(no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->parent[p] ]->nodeid);
//             }
//             fprintf(stderr, "\nCGraphnode child: ");
// 			for(int c=0; c<no2gnode->Get(i)->child.Count(); c++) {
//                 fprintf(stderr, "%d  ", no2gnode->Get(i)->child[c]); 
//                 fprintf(stderr, "%d  ", no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->child[c] ]->nodeid); 

//                 no2gnode_unispg[s][cgidx][i]->child.Add(no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->child[c] ]->nodeid);
//             }
//             fprintf(stderr, "\n");
//         }

//         // sink->parent.Add(node->nodeid); // add node to sink's parents

//         current_gidx[s] += 1;
//     } else {
//         fprintf(stderr, "************************************\n");
//         fprintf(stderr, "*********** Add more graph! ********\n");
//         fprintf(stderr, "************************************\n");
//         int next_unispg = true;
//         int graph_overlap = false;
//         while(next_unispg) {
//             int cgidx = current_gidx[s];
//             // Need to check overlapping & new graph creation

//             // for (int i = 0; i < no2gnode->Count(); i++) {
//             //     // no2gnode->Get(i);
//             //     fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
//             // }
//             // for (int i = 0; i < no2gnodeGp_unispg[s][cgidx].Count(); i++) {
//             //     fprintf(stderr, "&&&& This is the global graphnode: %d\n", no2gnodeGp_unispg[s][cgidx][i]->nodeid);
//             // }


//             uint lclg_start = no2gnode->Get(1)->start;
//             uint lclg_end = no2gnode->Get(no2gnode->Count()-2)->end;

//             uint unispg_start = no2gnode_unispg[s][cgidx][1]->start;
//             uint unispg_end = no2gnode_unispg[s][cgidx][ no2gnode_unispg[s][cgidx].Count()-2 ]->end;

//             fprintf(stderr, "$$$      cgidx: %d\n", cgidx);
//             fprintf(stderr, "$$$      no2gnode_unispg[s]: %d\n", sizeof(no2gnode_unispg[s]));
//             fprintf(stderr, "$$$ lclg_start: %u,  lclg_end: %u,  unispg_start: %u,  unispg_end: %u\n", lclg_start, lclg_end, no2gnode_unispg[s][cgidx][1]->start, no2gnode_unispg[s][cgidx][ no2gnode_unispg[s][cgidx].Count()-2 ]->end);
            
//             // unispg: -------
//             // lclg: ........
//             if (unispg_end < lclg_start) {
//                 // ----------   |(s).................(e)|
//                 fprintf(stderr,"\n  &&& Graph: ----------   |(s).................(e)| \n");
//                 next_unispg = true;
//                 graph_overlap = false;
//             } else if (unispg_start < lclg_start && unispg_end >= lclg_start && unispg_end <= lclg_end) {
//                 // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
//                 fprintf(stderr,"\n  &&& Graph: ----------|(s).................(e)|   or   -----|(s)-----............(e)| \n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start < lclg_start && unispg_end > lclg_end) {
//                 // -----|(s)------------(e)|--
//                 fprintf(stderr,"\n &&& Graph: -----|(s)------------(e)|-- \n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
//                 // |(s)----------.................(e)| 
//                 fprintf(stderr,"\n &&& Graph: |(s)----------............(e)|\n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start > lclg_start && unispg_end < lclg_end) {
//                 // |(s)........----------........(e)|
//                 fprintf(stderr,"\n &&& Graph: |(s)........----------........(e)| \n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start > lclg_start && unispg_end == lclg_end) {
//                 // |(s)............----------(e)|
//                 fprintf(stderr,"\n &&& Graph: |(s)............----------(e)| \n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
//                 // |(s)----------(e)|
//                 fprintf(stderr,"\n &&& Graph: |(s)----------(e)| \n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start <= lclg_end && unispg_end > lclg_end) {
//                 // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
//                 fprintf(stderr,"\n &&& Graph: (s)...............------(e)|-----    or   |(s).................(e)|---------- \n");
//                 next_unispg = true;
//                 graph_overlap = true;
//             } else if (unispg_start > lclg_end) {
//                 // The node is outside the current bundle => This node belongs to the next bundlenode
//                 // |(s).................(e)|   ----------
//                 fprintf(stderr,"\n &&& Graph: |(s).................(e)|   ---------- \n");
//                 next_unispg = false;
//                 graph_overlap = false;
//                 current_gidx[s] -= 1;
//             } else {
//                 fprintf(stderr,"\n &&& Unknown area!!!! \n");
//                 next_unispg = false;
//                 graph_overlap = false;
//             }
//             for (int i = 0; i < no2gnode->Count(); i++) {
//                 // no2gnode->Get(i);
//                 fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
//             }


//             /*****************************
//              * Create a new graph 
//              *****************************/
//             if (graph_overlap == true) {

//                 unsigned int nodeid = 0;
//                 // no2gnode_unispg[s][cgidx];

//                 GPVec<CGraphnodeUnispg>* new_no2gnode_unispg; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
// 			    new_no2gnode_unispg = new GPVec<CGraphnodeUnispg>[20000];


//                 int lclg_idx = 1;
//                 int unispg_idx = 1;
//                 bool more_comparison = true;
//                 bool lclg_is_end = false;
//                 bool unispg_is_end = false;
//                 while(!lclg_is_end || !unispg_is_end) {
//                     lclg_is_end = (lclg_idx == no2gnode->Count()-2);
//                     unispg_is_end = (unispg_idx == no2gnode_unispg[s][cgidx].Count()-2);

//                     CGraphnode* lclg_node =  no2gnode->Get(lclg_idx);
//                     fprintf(stderr, "&& lclg_node: %d (%d - %d)\t", lclg_node->nodeid, lclg_node->start, lclg_node->end);
//                     CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][cgidx].Get(unispg_idx);
//                     fprintf(stderr, "&& unispg_node: %d (%d - %d)\t", unispg_node->nodeid, unispg_node->start, unispg_node->end);


//                     bool lclg_move = false;
//                     bool unispg_move = false;
//                     // unispg: -------
//                     // lclg: |........|
//                     if (unispg_node->end < lclg_node->start) {
//                         fprintf(stderr,"\n  >>>>  Graph node: ----------   |(s).................(e)| \n\n");
// 	// CGraphnodeUnispg(int s=0,int e=0,unsigned int id=MAX_NODE, GVec<bool> is_passed_s_i=NULL, bool is_passed=false, GVec<float> cov_s_i=NULL, float cov=0, GVec<float> capacity_s_i=NULL, float capacity=0,float r=0):GSeg(s,e),nodeid(id),is_passed_s(),cov_s(),capacity_s(),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){

// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // To create a universal graph node, I should input a universal graphnode as well as a local graphnode.
// CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(unispg_node->start, unispg_node->end, nodeid, unispg_node->is_passed_s, false, unispg_node->cov_s, 0, unispg_node->capacity_s, 0, 0);
// CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(lclg_node->start, lclg_node->end, nodeid, NULL, true, NULL, lclg_node->cov, NULL, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);


//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                     } else if (unispg_node->start < lclg_node->start && unispg_node->end == lclg_node->start) {
//                         fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n\n");
// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(unispg_node->start, unispg_node->end, nodeid, false, 0, 0, 0);
// CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(lclg_node->start, lclg_node->end, nodeid, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);

//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                     } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->start && unispg_node->end < lclg_node->end) {
//                         fprintf(stderr,"\n  >>>>  Graph node: -----|(s)-----............(e)| \n\n");

// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, lclg_node->start, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(lclg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_3 = new CGraphnode(unispg_node->end, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                     } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->start && unispg_node->end == lclg_node->end) {
//                         fprintf(stderr,"\n  >>>>  Graph node: -----|(s)---------(e)| \n\n");

// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, lclg_node->start, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);

//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->end) {
//                         fprintf(stderr,"\n >>>>  Graph node: -----|(s)------------(e)|-- \n\n");


// 	// CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0):GSeg(s,e),
// 	// 		nodeid(id),cov(nodecov),capacity(cap),rate(r),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){}

// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, lclg_node->start, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_3 = new CGraphnode(lclg_node->end, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);


//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start == lclg_node->start && unispg_node->end < lclg_node->end) {
// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->end, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
//                         fprintf(stderr,"\n >>>>  Graph node: |(s)----------............(e)|\n\n");
//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                     } else if (unispg_node->start == lclg_node->start && unispg_node->end == lclg_node->end) {
// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
//                         fprintf(stderr,"\n >>>>  Graph node: |(s)----------(e)| \n\n");
//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start == lclg_node->start && unispg_node->end > lclg_node->end) {
//                         fprintf(stderr,"\n >>>>  Graph node: |(s)----------(e)|------ \n\n");
// // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->end, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start > lclg_node->start && unispg_node->end < lclg_node->end) {
//                         fprintf(stderr,"\n >>>>  Graph node: |(s)........----------........(e)| \n\n");
// // CGraphnode* node_1 = new CGraphnode(lclg_node->start, unispg_node->start, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_3 = new CGraphnode(unispg_node->end, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                     } else if (unispg_node->start > lclg_node->start && unispg_node->end == lclg_node->end) {
//                         fprintf(stderr,"\n >>>>  Graph node: |(s)............----------(e)| \n\n");
// // CGraphnode* node_1 = new CGraphnode(lclg_node->start, unispg_node->start, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);

//                         if (!unispg_is_end) {
//                             unispg_move = true;
//                         }
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start < lclg_node->start && unispg_node->start > lclg_node->end && unispg_node->end > lclg_node->end) {
//                         fprintf(stderr,"\n >>>>  Graph node: |(s)...............------(e)|----- \n\n");
// // CGraphnode* node_1 = new CGraphnode(lclg_node->start, unispg_node->start, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->start, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// // CGraphnode* node_3 = new CGraphnode(lclg_node->end, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start == lclg_node->end && unispg_node->end > lclg_node->end) {
//                         fprintf(stderr,"\n >>>>  Graph node: |(s).................(e)|---------- \n\n");
// // CGraphnode* node_1 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     } else if (unispg_node->start > lclg_node->end) {
//                         // The node is outside the current bundle => This node belongs to the next bundlenode
//                         fprintf(stderr,"\n >>>>  Graph node: |(s).................(e)|   ---------- \n\n");
// // CGraphnode* node_1 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
//                         if (!lclg_is_end) {
//                             lclg_move = true;
//                         }
//                     }
//                     if (lclg_move) {
//                         lclg_idx += 1;
//                     }
//                     if (unispg_move) {
//                         unispg_idx += 1;
//                     }
//                     if (!lclg_move && !unispg_move) {
//                         // Both lclg and unispg can not move forward.
//                         fprintf(stderr, ">>>> Local & global graph node cannot move\n");
//                         // fprintf(stderr, ">>>> Local & global graph node cannot move\n")
//                         // while(!lclg_is_end || !unispg_is_end) {
//                         //     lclg_is_end = (lclg_idx == no2gnode->Count()-2);
//                         // }
//                         if (!lclg_is_end) {
//                             lclg_idx += 1;
//                         }
//                         if (!unispg_is_end) {
//                             unispg_idx += 1;
//                         }
//                         if (lclg_is_end && unispg_is_end) {
//                             break;
//                         }
//                     }
//                 }
//             }
//             current_gidx[s] += 1;
//         }



//     }

//     // int g_idx = 0;
//     // for (int i=0; i<no2gnode->Count(); i++) {
//     //     no2gnodeGp_unispg[s][g_idx].Add(no2gnode->Get(i));
//     //     g_idx++;
//     // }


//     // gpSize[uni_splice_graph->get_s()] = uni_splice_graph->get_g_idx()+1;
//     // fprintf(stderr, "&&& uni_splice_graph->get_g_idx()+1: %d\n", uni_splice_graph->get_g_idx()+1);
//     // fprintf(stderr, "&&& gpSize[uni_splice_graph->get_s()]: %d\n", gpSize[uni_splice_graph->get_s()]);
//     // We need to reset graph_idx when (1)the new strands     
// }


// void Clear() {
//     // fprintf(stderr, "**** Start Clearing !!!! \n ");
//     for(int i=0;i<2;i++) {
//     // gpSize[i] = 0;
//     // graphnoGp[i].Clear();
//     // graphnoGp[i].setCapacity(8192);
//     // edgenoGp[i].Clear();
//     // edgenoGp[i].setCapacity(8192);
//     delete [] no2gnode_unispg[i];
//     no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[8192];
//     // no2gnodeGp_unispg[i]->setCapacity(8192);
//     };
// }

// // Need to be modifed!!!
// void PrintGraphGp() {
//     fprintf(stderr, "*********************************\n");
//     fprintf(stderr, "*********** PrintGraphGp ********\n");
//     fprintf(stderr, "*********************************\n");

//     { // DEBUG ONLY
//         printTime(stderr);
//         for(int s=0;s<2;s++) {
//             fprintf(stderr, "\n\tThere are %d stranded[%d] graphs\n", current_gidx[s]+1,int(2*s));
//             current_gidx[s];
//             for(int b=0;b<current_gidx[s]+1;b++) {                
//                 GStr pat;
//                 fprintf(stderr,"\t\tGraph[%d][%d] with %d nodes :",int(2*s),b,no2gnode_unispg[s][b].Count());
//                 for(int nd=1;nd<no2gnode_unispg[s][b].Count()-1;nd++) {
//                     fprintf(stderr, "&&& nd: %d\n", nd);
//                     fprintf(stderr, "&&& no2gnode_unispg[s][b][nd]->start: %u\n", no2gnode_unispg[s][b][nd]->start);
//                     fprintf(stderr, "&&& no2gnode_unispg[s][b][nd]->end: %u\n", no2gnode_unispg[s][b][nd]->end);
//                     fprintf(stderr," %d(%u-%u)",nd,no2gnode_unispg[s][b][nd]->start,no2gnode_unispg[s][b][nd]->end);
//                 }
//                 fprintf(stderr,"\n");
//             }
//         }

//     //     for(int s=0;s<2;s++) {
//     // 	fprintf(stderr, "There are %d stranded[%d] graphs\n",bno[s],int(2*s));
//     // 	// fprintf(stderr, "1 bundle[sno].Count(): %d\n", bno[s]);
//     // 	for(int b=0;b<bno[s];b++) {
//     // 		if(graphno[s][b]) {
//     // 			GStr pat;
//     // 			fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges with lastgpos=%d:",int(2*s),b,graphno[s][b],edgeno[s][b],lastgpos[s][b]);
//     // 			for(int nd=1;nd<graphno[s][b]-1;nd++)
//     // 				fprintf(stderr," %d(%d-%d)",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end);
//     // 			fprintf(stderr,"\n");
//     // 			print_pattern(tr2no[s][b],pat,graphno[s][b]);
//     // 		}
//     // 	}
//     // }
//     }
// }