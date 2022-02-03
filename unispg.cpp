#include "unispg.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

extern UnispgGp* unispg_gp;
extern GVec<int> current_gidx;

void UnispgGp::ProcessSample(GStr sample_name) {
    samples.Add(sample_name);
    for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
        int s=sno/2; // adjusted strand due to ignoring neutral strand
        current_gidx[s] = 0;
    }
}

void UnispgGp::AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode) {
    int sample_num = unispg_gp->samples.Count();
    if (fidx == 0) {
        fprintf(stderr, "*****************************\n");
        fprintf(stderr, "*********** AddGraph ********\n");
        fprintf(stderr, "*****************************\n");
        int cgidx = current_gidx[s];
        // This is the first unispg. Simply add copy it into UnispgGp
        for (int i=0; i<no2gnode->Count(); i++) {
            CGraphnode* node = new CGraphnode(no2gnode->Get(i));

			GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
			GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
			GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);

            CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(sample_num, node->start, node->end, node->nodeid, is_passed_s, cov_s, capacity_s, true, node->cov, node->capacity, node->rate);

            // no2gnodeGp_unispg[s][cgidx].Add(node);
            no2gnode_unispg[s][cgidx].Add(node_unispg);
        }

        // Linking parent and child
        for (int i=0; i<no2gnode->Count(); i++) {
            fprintf(stderr, "CGraphnode parent: ");
            for(int p=0; p<no2gnode->Get(i)->parent.Count(); p++) {
                fprintf(stderr, "%d  ", no2gnode->Get(i)->parent[p]);
                fprintf(stderr, "%d  ", no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->parent[p] ]->nodeid);

                no2gnode_unispg[s][cgidx][i]->parent.Add(no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->parent[p] ]);
            }
            fprintf(stderr, "\nCGraphnode child: ");
            for(int c=0; c<no2gnode->Get(i)->child.Count(); c++) {
                fprintf(stderr, "%d  ", no2gnode->Get(i)->child[c]); 
                fprintf(stderr, "%d  ", no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->child[c] ]->nodeid); 

                no2gnode_unispg[s][cgidx][i]->child.Add(no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->child[c] ]);
            }
            fprintf(stderr, "\n");
        }

        // sink->parent.Add(node->nodeid); // add node to sink's parents

        current_gidx[s] += 1;
    } else {
        fprintf(stderr, "************************************\n");
        fprintf(stderr, "*********** Add more graph! ********\n");
        fprintf(stderr, "************************************\n");
        int next_unispg = true;
        int process_unispg = false;
        int cgidx = current_gidx[s];
        fprintf(stderr, "More than 1 graph: %d\n", cgidx);

        while(next_unispg) {
            int cgidx = current_gidx[s];
            // Need to check overlapping & new graph creation

            // for (int i = 0; i < no2gnode->Count(); i++) {
            //     // no2gnode->Get(i);
            //     fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
            // }
            // for (int i = 0; i < no2gnodeGp_unispg[s][cgidx].Count(); i++) {
            //     fprintf(stderr, "&&&& This is the global graphnode: %d\n", no2gnodeGp_unispg[s][cgidx][i]->nodeid);
            // }


            uint lclg_start = no2gnode->Get(1)->start;
            uint lclg_end = no2gnode->Get(no2gnode->Count()-2)->end;

            uint unispg_start = no2gnode_unispg[s][cgidx][1]->start;
            uint unispg_end = no2gnode_unispg[s][cgidx][ no2gnode_unispg[s][cgidx].Count()-2 ]->end;

            fprintf(stderr, "$$$      cgidx: %d\n", cgidx);
            fprintf(stderr, "$$$      no2gnode_unispg[s]: %d\n", sizeof(no2gnode_unispg[s]));
            fprintf(stderr, "$$$ lclg_start: %u,  lclg_end: %u,  unispg_start: %u,  unispg_end: %u\n", lclg_start, lclg_end, no2gnode_unispg[s][cgidx][1]->start, no2gnode_unispg[s][cgidx][ no2gnode_unispg[s][cgidx].Count()-2 ]->end);
            
            // unispg: -------
            // lclg: ........
            if (unispg_end < lclg_start) {
                // ----------   |(s).................(e)|
                fprintf(stderr,"\n  &&& Graph: ----------   |(s).................(e)| \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start < lclg_start && unispg_end >= lclg_start && unispg_end <= lclg_end) {
                // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
                fprintf(stderr,"\n  &&& Graph: ----------|(s).................(e)|   or   -----|(s)-----............(e)| \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start < lclg_start && unispg_end > lclg_end) {
                // -----|(s)------------(e)|--
                fprintf(stderr,"\n &&& Graph: -----|(s)------------(e)|-- \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
                // |(s)----------.................(e)| 
                fprintf(stderr,"\n &&& Graph: |(s)----------............(e)|\n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start > lclg_start && unispg_end < lclg_end) {
                // |(s)........----------........(e)|
                fprintf(stderr,"\n &&& Graph: |(s)........----------........(e)| \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start > lclg_start && unispg_end == lclg_end) {
                // |(s)............----------(e)|
                fprintf(stderr,"\n &&& Graph: |(s)............----------(e)| \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
                // |(s)----------(e)|
                fprintf(stderr,"\n &&& Graph: |(s)----------(e)| \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start <= lclg_end && unispg_end > lclg_end) {
                // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
                fprintf(stderr,"\n &&& Graph: (s)...............------(e)|-----    or   |(s).................(e)|---------- \n");
                next_unispg = true;
                process_unispg = true;
            } else if (unispg_start > lclg_end) {
                // The node is outside the current bundle => This node belongs to the next bundlenode
                // |(s).................(e)|   ----------
                fprintf(stderr,"\n &&& Graph: |(s).................(e)|   ---------- \n");
                next_unispg = false;
                process_unispg = false;
                current_gidx[s] -= 1;
            } else {
                fprintf(stderr,"\n &&& Unknown area!!!! \n");
                next_unispg = false;
                process_unispg = false;
            }
            for (int i = 0; i < no2gnode->Count(); i++) {
                // no2gnode->Get(i);
                fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
            }


            /*****************************
             * Create a new graph 
             *****************************/
            if (process_unispg) {

                unsigned int unispg_nodeid = 1;
                // no2gnode_unispg[s][cgidx];
                GPVec<CGraphnodeUnispg>* new_no2gnode_unispg; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
                new_no2gnode_unispg = new GPVec<CGraphnodeUnispg>[20000];


                int lclg_idx = 1;
                int unispg_idx = 1;

                int prev_lclg_idx = 1;
                int prev_unispg_idx = 1;

                bool more_comparison = true;
                bool lclg_is_end = false;
                bool unispg_is_end = false;

	            GHashMap<int, GVec<CGraphnodeUnispg*>* > lclg_2_unispg(false); //hash of pointers
                // lclg_2_unispg

                GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, is_passed_s, cov_s, capacity_s, false, 0, 0, 0);

                new_no2gnode_unispg[cgidx].Add(source);
        

// need to relabel the nodeid
                while(!lclg_is_end || !unispg_is_end) {
                    lclg_is_end = (lclg_idx == no2gnode->Count()-2);
                    unispg_is_end = (unispg_idx == no2gnode_unispg[s][cgidx].Count()-2);

                    CGraphnode* lclg_node =  no2gnode->Get(lclg_idx);
                    fprintf(stderr, "&& lclg_node: %d (%d - %d)\t", lclg_node->nodeid, lclg_node->start, lclg_node->end);
                    CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][cgidx].Get(unispg_idx);
                    fprintf(stderr, "&& unispg_node: %d (%d - %d)\t", unispg_node->nodeid, unispg_node->start, unispg_node->end);


                    bool lclg_move = false;
                    bool unispg_move = false;
                    // unispg: -------
                    // lclg: |........|
                    if (unispg_node->end < lclg_node->start) {
                        fprintf(stderr,"\n  >>>>  Graph node: ----------   |(s).................(e)| \n\n");
    // CGraphnodeUnispg(int s=0,int e=0,unsigned int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, bool is_passed=false, GVec<float>* cov_s_i=NULL, float cov=0, GVec<float>* capacity_s_i=NULL, float capacity=0,float r=0):GSeg(s,e),nodeid(id),is_passed_s(),cov_s(),capacity_s(),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){

// CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// To create a universal graph node, I should input a universal graphnode as well as a local graphnode.
// fprintf(stderr, "CGraphnode parent: ");
// for(int p=0; p<unispg_node->parent.Count(); p++) {
//     fprintf(stderr, "%d  ", unispg_node->parent[p]->nodeid);
//     fprintf(stderr, "%d  ", unispg_gp->no2gnode_unispg[s][cgidx][ unispg_node->parent[p]->nodeid ]->nodeid);

//     new_no2gnode_unispg[cgidx][unispg_nodeid]->parent.Add(unispg_node->parent[p]);
// }
// fprintf(stderr, "\nCGraphnode child: ");
// for(int c=0; c<unispg_node->child.Count(); c++) {
//     fprintf(stderr, "%d  ", unispg_node->child[c]->nodeid);
//     fprintf(stderr, "%d  ", unispg_gp->no2gnode_unispg[s][cgidx][ unispg_node->child[c]->nodeid ]->nodeid);

//     new_no2gnode_unispg[cgidx][unispg_nodeid]->child.Add(unispg_node->child[c]);
// }
// fprintf(stderr, "\n");

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;


GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);

// for (int i = 0; i < new_is_passed_s->Count(); i++) {
//     fprintf(stderr, "Outside: %d\n", new_is_passed_s->Get(i));
// }
// for (int i = 0; i < new_cov_s->Count(); i++) {
//     fprintf(stderr, "Outside: %f\n", new_cov_s->Get(i));
// }
// for (int i = 0; i < new_capacity_s->Count(); i++) {
//     fprintf(stderr, "Outside: %f\n", new_capacity_s->Get(i));
// }
CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
// CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;


GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);
fprintf(stderr, "Checker!!!: %d\n", lclg_2_unispg[lclg_node->nodeid]->Get(0)->nodeid);

// fprintf(stderr, "CGraphnode parent: ");
// for(int p=0; p<lclg_node->parent.Count(); p++) {
//     fprintf(stderr, "%d  ", lclg_node->parent[p]);
//     // fprintf(stderr, "%d  ", unispg_gp->no2gnode_unispg[s][cgidx][ unispg_node->parent[p]->nodeid ]->nodeid);

//     // new_no2gnode_unispg[cgidx][unispg_nodeid]->parent.Add(unispg_node->parent[p]);
// }
// fprintf(stderr, "\nCGraphnode child: ");
// for(int c=0; c<lclg_node->child.Count(); c++) {
//     fprintf(stderr, "%d  ", lclg_node->child[c]);
//     // fprintf(stderr, "%d  ", unispg_gp->no2gnode_unispg[s][cgidx][ unispg_node->child[c]->nodeid ]->nodeid);

//     // new_no2gnode_unispg[cgidx][unispg_nodeid]->child.Add(unispg_node->child[c]);
// }
// fprintf(stderr, "\n");

                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                    } else if (unispg_node->start < lclg_node->start && unispg_node->end == lclg_node->start) {
                        fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n\n");
// CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
// CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;


GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                    } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->start && unispg_node->end < lclg_node->end) {
                        fprintf(stderr,"\n  >>>>  Graph node: -----|(s)-----............(e)| \n\n");

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->start, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_3 = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_3);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_2);
nodelist->Add(node_3);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


// CGraphnode* node_3 = new CGraphnode(unispg_node->end, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                    } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->start && unispg_node->end == lclg_node->end) {
                        fprintf(stderr,"\n  >>>>  Graph node: -----|(s)---------(e)| \n\n");

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->start, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);

                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: -----|(s)------------(e)|-- \n\n");

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->start, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0,  0, 0);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

CGraphnodeUnispg* node_3 = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0,  0, 0);
new_no2gnode_unispg[cgidx].Add(node_3);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start == lclg_node->start && unispg_node->end < lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s)----------............(e)|\n\n");


CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                    } else if (unispg_node->start == lclg_node->start && unispg_node->end == lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s)----------(e)| \n\n");

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start == lclg_node->start && unispg_node->end > lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s)----------(e)|------ \n\n");

CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);

                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start > lclg_node->start && unispg_node->end < lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s)........----------........(e)| \n\n");

GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->start, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<bool>* new_is_passed_s_2 = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s_2 = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s_2 = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_3 = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, unispg_nodeid, new_is_passed_s_2, new_cov_s_2, new_capacity_s_2, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_3);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
nodelist->Add(node_2);
nodelist->Add(node_3);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);

                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                    } else if (unispg_node->start > lclg_node->start && unispg_node->end == lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s)............----------(e)| \n\n");


GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->start, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


                        if (!unispg_is_end) {
                            unispg_move = true;
                        }
                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start < lclg_node->start && unispg_node->start > lclg_node->end && unispg_node->end > lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s)...............------(e)|----- \n\n");

GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->start, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

CGraphnodeUnispg* node_3 = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_3);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
nodelist->Add(node_2);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);


                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start == lclg_node->end && unispg_node->end > lclg_node->end) {
                        fprintf(stderr,"\n >>>>  Graph node: |(s).................(e)|---------- \n\n");

GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);

                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    } else if (unispg_node->start > lclg_node->end) {
                        // The node is outside the current bundle => This node belongs to the next bundlenode
                        fprintf(stderr,"\n >>>>  Graph node: |(s).................(e)|   ---------- \n\n");

GVec<bool>* new_is_passed_s = new GVec<bool>(sample_num-1, false);
GVec<float>* new_cov_s = new GVec<float>(sample_num-1, 0.0f);
GVec<float>* new_capacity_s = new GVec<float>(sample_num-1, 0.0f);
CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, unispg_nodeid, new_is_passed_s, new_cov_s, new_capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
new_no2gnode_unispg[cgidx].Add(node_1);
unispg_nodeid += 1;

CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
new_no2gnode_unispg[cgidx].Add(node_2);
unispg_nodeid += 1;

GVec<CGraphnodeUnispg*>* nodelist = new GVec<CGraphnodeUnispg*>[1];
nodelist->Add(node_1);
lclg_2_unispg.Add(lclg_node->nodeid, nodelist);

                        if (!lclg_is_end) {
                            lclg_move = true;
                        }
                    }
                    if (lclg_move) {
                        lclg_idx += 1;
                    }
                    if (unispg_move) {
                        unispg_idx += 1;
                    }
                    if (!lclg_move && !unispg_move) {
                        // Both lclg and unispg can not move forward.
                        fprintf(stderr, ">>>> Local & global graph node cannot move\n");
                        // fprintf(stderr, ">>>> Local & global graph node cannot move\n")
                        // while(!lclg_is_end || !unispg_is_end) {
                        //     lclg_is_end = (lclg_idx == no2gnode->Count()-2);
                        // }
                        if (!lclg_is_end) {
                            lclg_idx += 1;
                        }
                        if (!unispg_is_end) {
                            unispg_idx += 1;
                        }
                        if (lclg_is_end && unispg_is_end) {
                            break;
                        }
                    }
                    prev_lclg_idx = lclg_idx;
                    prev_unispg_idx = unispg_idx;
                }
            }
            current_gidx[s] += 1;
        }
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
    // 				fprintf(stderr," %d(%d-%d)",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end);
    // 			fprintf(stderr,"\n");
    // 			print_pattern(tr2no[s][b],pat,graphno[s][b]);
    // 		}
    // 	}
    // }
    }
}

GPVec<CGraphnodeUnispg>** UnispgGp::get_no2gnodeGp () {
    return no2gnode_unispg;
}