#include "unispg.h"
#include "GBitVec.h"
#include <float.h>

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

void UnispgGp::AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode_base, int lclg_idx, int lclg_limit) {
    // Next lclg: no2gnode+1 !!    
    int sample_num = unispg_gp->samples.Count();

    GPVec<CGraphnode>* no2gnode;
    int unispg_idx = current_gidx[s];
    if (fidx == 0) {
        fprintf(stderr, "\n*****************************\n");
        fprintf(stderr, "*********** AddGraph ********\n");
        fprintf(stderr, "*****************************\n");
        no2gnode = no2gnode_base + lclg_idx;
        // This is the first unispg. Simply add copy it into UnispgGp
        for (int i=0; i<no2gnode->Count(); i++) {
            CGraphnode* node = new CGraphnode(no2gnode->Get(i));

			GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
			GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
			GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);

            CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(sample_num, node->start, node->end, node->nodeid, is_passed_s, cov_s, capacity_s, true, node->cov, node->capacity, node->rate);

            // no2gnodeGp_unispg[s][unispg_idx].Add(node);
            no2gnode_unispg[s][unispg_idx].Add(node_unispg);
        }

        // Linking parent and child
        for (int i=0; i<no2gnode->Count(); i++) {
            fprintf(stderr, "CGraphnode parent: ");
            for(int p=0; p<no2gnode->Get(i)->parent.Count(); p++) {
                fprintf(stderr, "%d  ", no2gnode->Get(i)->parent[p]);
                fprintf(stderr, "%d  ", no2gnode_unispg[s][unispg_idx][ no2gnode->Get(i)->parent[p] ]->nodeid);

                no2gnode_unispg[s][unispg_idx][i]->parent.Add(no2gnode_unispg[s][unispg_idx][ no2gnode->Get(i)->parent[p] ]);
            }
            fprintf(stderr, "\nCGraphnode child: ");
            for(int c=0; c<no2gnode->Get(i)->child.Count(); c++) {
                fprintf(stderr, "%d  ", no2gnode->Get(i)->child[c]); 
                fprintf(stderr, "%d  ", no2gnode_unispg[s][unispg_idx][ no2gnode->Get(i)->child[c] ]->nodeid); 

                no2gnode_unispg[s][unispg_idx][i]->child.Add(no2gnode_unispg[s][unispg_idx][ no2gnode->Get(i)->child[c] ]);
            }
            fprintf(stderr, "\n");
        }

        // sink->parent.Add(node->nodeid); // add node to sink's parents

        current_gidx[s] += 1;
    } else {
        fprintf(stderr, "************************************\n");
        fprintf(stderr, "*********** Add more graph! ********\n");
        fprintf(stderr, "************************************\n");

        bool more_unispg = true;
        bool more_lclg = true;
        int next_unispg = true;
        int process_unispg = false;
        // int lclg_idx = 0;
        fprintf(stderr, "More than 1 graph: %d\n", unispg_idx);

        // Here, I need to find out how many unispg & lclg should be merged into 1 new unispg.

        GPVec<CGraphnodeUnispg>* unispg_merge = new GPVec<CGraphnodeUnispg>[200];
        GPVec<CGraphnode>* lclg_merge = new GPVec<CGraphnode>[200];
        int unispg_merge_idx = 0;
        int lclg_merge_idx = 0;

        while(more_unispg || more_lclg) {
            fprintf(stderr, "unispg_merge_idx: %d, lclg_idx: %d, lclg_limit: %d\n", unispg_merge_idx, lclg_merge_idx, lclg_limit);
            unispg_idx = current_gidx[s];
            no2gnode = no2gnode_base + lclg_idx;

            fprintf(stderr, "no2gnode: %p\n", no2gnode);
            fprintf(stderr, "count: %d \n", no2gnode->Count());
            uint lclg_start = no2gnode->Get(1)->start;
            uint lclg_end = no2gnode->Get(no2gnode->Count()-2)->end;
            fprintf(stderr, "boundary: %d - %d \n", lclg_start, lclg_end);

            uint unispg_start = no2gnode_unispg[s][unispg_idx][1]->start;
            uint unispg_end = no2gnode_unispg[s][unispg_idx][ no2gnode_unispg[s][unispg_idx].Count()-2 ]->end;

            fprintf(stderr, "$$$      unispg_idx: %d\n", unispg_idx);
            fprintf(stderr, "$$$      no2gnode_unispg[s]: %d\n", sizeof(no2gnode_unispg[s]));
            fprintf(stderr, "$$$ lclg_start: %u,  lclg_end: %u,  unispg_start: %u,  unispg_end: %u\n", lclg_start, lclg_end, no2gnode_unispg[s][unispg_idx][1]->start, no2gnode_unispg[s][unispg_idx][ no2gnode_unispg[s][unispg_idx].Count()-2 ]->end);
            

            for (int i = 0; i < no2gnode->Count(); i++) {
                // no2gnode->Get(i);
                fprintf(stderr, "Before ~~ &&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
            }

            // unispg: -------
            // lclg: ........
            if (unispg_end < lclg_start) {
                // ----------   |(s).................(e)|
                fprintf(stderr,"\n  &&& Graph: ----------   |(s).................(e)|\n");
                // This is the end of the lclg & unispg comparison. Only need to use the unispg. 
                more_unispg = true;
                more_lclg = false;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                process_unispg = true; // Process when it is the end
            } else if (unispg_end == lclg_start) {
                // ----------|(s).................(e)|
                fprintf(stderr,"\n  &&& Graph: ----------|(s).................(e)| \n");
                more_unispg = true;
                more_lclg = false;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                process_unispg = false;
            } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end < lclg_end) {
                // -----|(s)-----............(e)|
                fprintf(stderr,"\n  &&& Graph: -----|(s)-----............(e)| \n");
                more_unispg = true;
                more_lclg = false;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                process_unispg = false;
            } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end == lclg_end) {
                // -----|(s)--------------(e)|
                fprintf(stderr,"\n  &&& Graph: -----|(s)--------------(e)| \n");
                //  This is the end of the lclg & unispg comparison. Set false first, but I need to check whether the next unispg and lclg overlap with the current graph.
                more_unispg = true;
                more_lclg = true;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = true; // Process when it is the end
            } else if (unispg_start < lclg_start && unispg_end > lclg_start && unispg_end > lclg_end) {
                // -----|(s)------------(e)|--
                fprintf(stderr,"\n &&& Graph: -----|(s)------------(e)|-- \n");
                more_unispg = false;
                more_lclg = true;
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = false;
            } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
                // |(s)----------.................(e)| 
                fprintf(stderr,"\n &&& Graph: |(s)----------............(e)|\n");
                more_unispg = true;
                more_lclg = false;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                process_unispg = false;
            } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
                // |(s)----------(e)|
                fprintf(stderr,"\n &&& Graph: |(s)----------(e)| \n");
                // This is the end of the lclg & unispg comparison. 
                more_unispg = true;
                more_lclg = true;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = true; // Process when it is the end
            } else if (unispg_start == lclg_start && unispg_end > lclg_end) {
                // |(s)----------(e)|-----
                fprintf(stderr,"\n &&& Graph: |(s)----------(e)|----- \n");
                more_unispg = false;
                more_lclg = true;
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = false;
            } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end < lclg_end) {
                // |(s)........----------........(e)|
                fprintf(stderr,"\n &&& Graph: |(s)........----------........(e)| \n");
                more_unispg = true;
                more_lclg = false;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                process_unispg = false;
            } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end == lclg_end) {
                // |(s)............----------(e)|
                fprintf(stderr,"\n &&& Graph: |(s)............----------(e)| \n");
                // This is the end of the lclg & unispg comparison. 
                more_unispg = true;
                more_lclg = true;
                unispg_merge[unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = true; // Process when it is the end
            } else if (unispg_start > lclg_start && unispg_start < lclg_end && unispg_end > lclg_end) {
                // |(s)...............------(e)|-----
                fprintf(stderr,"\n &&& Graph: |(s)...............------(e)|-----  \n");
                more_unispg = false;
                more_lclg = true;
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = false;
            } else if (unispg_start == lclg_end) {
                // The node is outside the current bundle => This node belongs to the next bundlenode
                // |(s).................(e)|----------
                fprintf(stderr,"\n &&& Graph: |(s).................(e)|----------  \n");
                more_unispg = false;
                more_lclg = true;
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = false;
            } else if (unispg_start > lclg_end) {
                // The node is outside the current bundle => This node belongs to the next bundlenode
                // |(s).................(e)|   ----------
                fprintf(stderr,"\n &&& Graph: |(s).................(e)|   ---------- \n");
                // This is the end of the lclg & unispg comparison. Only need to use the lclg.
                more_unispg = false;
                more_lclg = true;
                // [unispg_merge_idx] = no2gnode_unispg[s][unispg_idx];
                lclg_merge[lclg_merge_idx] = no2gnode_base[lclg_idx];
                process_unispg = true; // Process when it is the end
            } else {
                fprintf(stderr,"\n &&& Unknown area!!!! \n");
                next_unispg = false;
                process_unispg = false;
            }
            for (int i = 0; i < no2gnode->Count(); i++) {
                // no2gnode->Get(i);
                fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
            }

            // Check lclg bounary
            if (lclg_idx == lclg_limit-1) {
                process_unispg = true;
            }









            /*****************************
             * Create a new graph 
             *****************************/
            if (process_unispg) {
                fprintf(stderr, "Inside `process_unispg`!!!\n");
                for (int i = 0; i <= unispg_merge_idx; i++) {
                    fprintf(stderr, "unispg_merge -> i: %d\n", i);
                    for (int j = 0; j < unispg_merge[i].Count(); j++) {
                        fprintf(stderr, "** unispg_merge[i][j]->start: %d\n", unispg_merge[i][j]->start);
                        fprintf(stderr, "** unispg_merge[i][j]->end: %d\n", unispg_merge[i][j]->end);
                    }
                }
                for (int i = 0; i <= lclg_merge_idx; i++) {
                    fprintf(stderr, "lclg_merge_idx -> i: %d\n", i);
                    for (int j = 0; j < lclg_merge[i].Count(); j++) {
                        fprintf(stderr, "** lclg_merge[i][j]->start: %d\n", lclg_merge[i][j]->start);
                        fprintf(stderr, "** lclg_merge[i][j]->end: %d\n", lclg_merge[i][j]->end);
                    }
                }
                fprintf(stderr, "\n\n\n");

                unispg_merge->Clear();
                unispg_merge_idx = 0;
                lclg_merge->Clear();
                lclg_merge_idx = 0;
                // for (int i = 0; i < lclg_merge_idx; i++) {
                //     lclg_merge
                // }
                // GPVec<CGraphnode>* unispg_merge = new GPVec<CGraphnode>[200];
                // GPVec<CGraphnode>* lclg_merge = new GPVec<CGraphnode>[200];
//                 // Creating boundary list!
//                 GVec<uint> boundaries;
//                 GVec<CGraphBoundaryType> boundaries_types;

//                 unsigned int new_unispg_nodeid = 1;
//                 // no2gnode_unispg[s][unispg_idx];
//                 GPVec<CGraphnodeUnispg>* new_no2gnode_unispg; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
//                 new_no2gnode_unispg = new GPVec<CGraphnodeUnispg>[20000];

//                 int lclg_idx = 1;
//                 int unispg_idx = 1;

//                 bool more_comparison = true;
//                 bool lclg_is_end = false;
//                 bool unispg_is_end = false;

//                 GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                 GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                 GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                 CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, is_passed_s, cov_s, capacity_s, false, 0, 0, 0);

//                 new_no2gnode_unispg[unispg_idx].Add(source);

//                 // Iterating boundaries.
//                 // int cur_boundary_start = 0;
//                 // int cur_boundary_end = 0;
//                 // int prev_boundary_start = 0;
//                 int prev_boundary = 0;

//                 CGraphnode* lclg_node =  no2gnode->Get(1);
//                 CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][unispg_idx].Get(1);
//                 prev_boundary = (unispg_node->start < lclg_node->start) ? unispg_node->start:lclg_node->start;
//                 while(!lclg_is_end || !unispg_is_end) {
//                     lclg_is_end = (lclg_idx == no2gnode->Count()-2);
//                     unispg_is_end = (unispg_idx == no2gnode_unispg[s][unispg_idx].Count()-2);

//                     CGraphnode* lclg_node =  no2gnode->Get(lclg_idx);
//                     // fprintf(stderr, "&& lclg_node: %d (%d - %d)\t", lclg_node->nodeid, lclg_node->start, lclg_node->end);
//                     CGraphnodeUnispg* unispg_node = no2gnode_unispg[s][unispg_idx].Get(unispg_idx);
//                     // fprintf(stderr, "&& unispg_node: %d (%d - %d)\t", unispg_node->nodeid, unispg_node->start, unispg_node->end);
//                     CGraphnodeUnispg* node;

//                     fprintf(stderr, "****** >> prev_boundary: %d \n", prev_boundary);
//                     fprintf(stderr, "****** >> lclg_idx: %d,  unispg_idx: %d\n", lclg_idx, unispg_idx);
//                     fprintf(stderr, "****** >> lclg_node->start: %d,  lclg_node->end: %d\n", lclg_node->start, lclg_node->end);
//                     fprintf(stderr, "****** >> unispg_node->start: %d,  unispg_node->end: %d\n", unispg_node->start, unispg_node->end);

//                     bool lclg_move = false;
//                     bool unispg_move = false;

//                     if (unispg_node->start < lclg_node->start) {
//                         AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
//                         uint node_start_pos = 0;
//                         if (prev_boundary < unispg_node->start) {
//                             node_start_pos = unispg_node->start;
//                         } else if (prev_boundary == unispg_node->start) {
//                             node_start_pos = unispg_node->start;
//                         } else if (prev_boundary > unispg_node->start) {
//                             node_start_pos = prev_boundary;
//                         }
//                         if (unispg_node->end < lclg_node->start) {
//                             fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
//                             AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             prev_boundary = unispg_node->end;
//                             MoveUnispg(unispg_is_end, unispg_move);
//                             // if (unispg_is_end) {
//                             //     AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
//                             //     AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                             // }
//                         } else if (unispg_node->end == lclg_node->start) {
//                             fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
//                             AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_S);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             prev_boundary = unispg_node->end;
//                             MoveUnispg(unispg_is_end, unispg_move);
//                         } else if (lclg_node->start < unispg_node->end) {
//                             AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->start, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             if (unispg_node->end < lclg_node->end) {
//                                 fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
//                                 AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;

//                                 prev_boundary = unispg_node->end;
//                                 MoveUnispg(unispg_is_end, unispg_move);
//                             } else if (unispg_node->end == lclg_node->end) {
//                                 fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----------(e)|\n");
//                                 AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, lclg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;


//                                 prev_boundary = unispg_node->end;
//                                 MoveUnispg(unispg_is_end, unispg_move);
//                                 MoveLclg(lclg_is_end, lclg_move);
//                             } else if (lclg_node->end < unispg_node->end) {
//                                 fprintf(stderr,"\n  >>>>  Graph node: -----|(s)----------(e)|-----\n");
//                                 AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;

//                                 prev_boundary = lclg_node->end;
//                                 MoveLclg(lclg_is_end, lclg_move);
//                             }
//                         }
//                     } else if (unispg_node->start == lclg_node->start) {
//                         AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S_LCLG_S);
//                         uint node_start_pos = 0;
//                         if (prev_boundary < unispg_node->start) {
//                             node_start_pos = unispg_node->start;
//                         } else if (prev_boundary == unispg_node->start) {
//                             node_start_pos = unispg_node->start;
//                         } else if (prev_boundary > unispg_node->start) {
//                             node_start_pos = prev_boundary;
//                         }


//                         if (unispg_node->end < lclg_node->end) {
//                             fprintf(stderr,"\n  >>>>  Graph node: |(s)----------.......(e)|\n");
//                             AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             prev_boundary = unispg_node->end;
//                             MoveUnispg(unispg_is_end, unispg_move);
//                         } else if (unispg_node->end == lclg_node->end) {
//                             fprintf(stderr,"\n  >>>>  Graph node: |(s)------------(e)|\n");
//                             AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             prev_boundary = unispg_node->end;
//                             MoveUnispg(unispg_is_end, unispg_move);
//                             MoveLclg(lclg_is_end, lclg_move);
//                         } else if (lclg_node->end < unispg_node->end) {
//                             fprintf(stderr,"\n  >>>>  Graph node: |(s)---------(e)|---\n");
//                             AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             prev_boundary = lclg_node->end;
//                             MoveLclg(lclg_is_end, lclg_move);
//                         }
//                     } else if (lclg_node->start < unispg_node->start) {
//                         AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
//                         uint node_start_pos = 0;
//                         if (prev_boundary < lclg_node->start) {
//                             node_start_pos = lclg_node->start;
//                         } else if (prev_boundary == lclg_node->start) {
//                             node_start_pos = lclg_node->start;
//                         } else if (prev_boundary > lclg_node->start) {
//                             node_start_pos = prev_boundary;
//                         }
//                         if (unispg_node->start < lclg_node->end) {
//                             AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
//                             GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                             GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                             GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                             // Add new node
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->start, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;

//                             if (unispg_node->end < lclg_node->end) {
//                                 fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------....(e)|\n");
//                                 AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;

//                                 prev_boundary = unispg_node->end;
//                                 MoveUnispg(unispg_is_end, unispg_move);
//                             } else if (unispg_node->end == lclg_node->end) {
//                                 fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|\n");
//                                 AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E_LCLG_E);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;

//                                 prev_boundary = unispg_node->end;
//                                 MoveUnispg(unispg_is_end, unispg_move);
//                                 MoveLclg(lclg_is_end, lclg_move);
//                             } else if (lclg_node->end < unispg_node->end) {
//                                 fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
//                                 AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, unispg_node->start, lclg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;
                                
//                                 prev_boundary = lclg_node->end;
//                                 MoveLclg(lclg_is_end, lclg_move);
//                             }
//                         } else if (lclg_node->end == unispg_node->start) {
//                             fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|----------\n");
//                             AddBoundary(boundaries, lclg_node->end, boundaries_types, UNISPG_S_LCLG_E);
//                             uint node_start_pos = 0;
//                             bool create_node = true;
//                             if (prev_boundary < lclg_node->start) {
//                                 node_start_pos = lclg_node->start;
//                             } else if (prev_boundary == lclg_node->start) {
//                                 node_start_pos = lclg_node->start;
//                             } else if (prev_boundary > lclg_node->start) {
//                                 if (prev_boundary < lclg_node->end) {
//                                     node_start_pos = prev_boundary;
//                                 } else if (prev_boundary == lclg_node->end) {
//                                     // Do not need to create node
//                                     create_node = false;
//                                 } else if (prev_boundary > lclg_node->end) {
//                                     // This is an impossible case. Insane check
//                                 }
//                             }
//                             if (create_node) {
//                                 GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                 GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                 GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;
//                             }

//                             prev_boundary = lclg_node->end;
//                             MoveLclg(lclg_is_end, lclg_move);
//                         } else if (lclg_node->end < unispg_node->start) {
//                             fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n");
//                             AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                             uint node_start_pos = 0;
//                             bool create_node = true;
//                             if (prev_boundary < lclg_node->start) {
//                                 node_start_pos = lclg_node->start;
//                             } else if (prev_boundary == lclg_node->start) {
//                                 node_start_pos = lclg_node->start;
//                             } else if (prev_boundary > lclg_node->start) {
//                                 if (prev_boundary < lclg_node->end) {
//                                     node_start_pos = prev_boundary;
//                                 } else if (prev_boundary == lclg_node->end) {
//                                     // Do not need to create node
//                                     create_node = false;
//                                 } else if (prev_boundary > lclg_node->end) {
//                                     // This is an impossible case. Insane check
//                                 }
//                             }
//                             if (create_node) {
//                                 GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                 GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                 GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                 // Add new node
//                                 node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                 new_no2gnode_unispg[unispg_idx].Add(node);
//                                 new_unispg_nodeid += 1;
//                             }

//                             prev_boundary = lclg_node->end;
//                             MoveLclg(lclg_is_end, lclg_move);
//                             // if (lclg_is_end) {
//                             //     AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
//                             //     AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                             // }
//                         }
//                     }

//                     if (lclg_move) {
//                         fprintf(stderr, "lclg_idx Move ! %d\n", lclg_idx);
//                         lclg_idx += 1;
//                     }
//                     if (unispg_move) {
//                         fprintf(stderr, "unispg_idx Move ! %d\n", unispg_idx);
//                         unispg_idx += 1;
//                     }
//                     if (!lclg_move && !unispg_move) {
//                         // Both lclg and unispg can not move forward.
//                         fprintf(stderr, ">>>> Local & global graph node cannot move\n");
//                         // fprintf(stderr, ">>>> Local & global graph node cannot move\n")
//                         // while(!lclg_is_end || !unispg_is_end) {
//                         //     lclg_is_end = (lclg_idx == no2gnode->Count()-2);
//                         // }

//                         if (lclg_is_end && !unispg_is_end) {
//                             AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
//                             AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                             uint node_start_pos = 0;
//                             if (prev_boundary < lclg_node->end) {
//                                 // it is impossible. Just the insane check
//                                 GError("Wrong boundaries. Check the code!!");
//                             } else if (prev_boundary == lclg_node->end) {
//                                 if (unispg_node->start < prev_boundary) {
// fprintf(stderr,"\n  >>>>  Graph node: -----|(s)----------(e)|-----\n");
// fprintf(stderr,"\n  >>>>  Graph node: |(s)---------(e)|---\n");
// fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
//                                     node_start_pos = prev_boundary;
//                                 } else if (unispg_node->start >= prev_boundary) {
// fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|----------\n");
// fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n");
//                                     node_start_pos = unispg_node->start;
//                                 }
//                             } else if (prev_boundary > lclg_node->end) {
//                                 // This situation has happend at least one. Insert the new unispg.
//                                 node_start_pos = unispg_node->start;
//                             }
//                             // Add new node
//                             fprintf(stderr, "^^^^ node_start_pos: %d,  unispg_node->end: %d \n", node_start_pos, unispg_node->end);
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;
//                             // !!!!!! Do not update the prev_boundary!!!
//                             // prev_boundary = unispg_node->end;
//                             // Move to next node of the unispg
//                             unispg_idx += 1;
//                         } else if (!lclg_is_end && unispg_is_end) {
//                             AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
//                             AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);

//                             uint node_start_pos = 0;
//                             if (prev_boundary < unispg_node->end) {
//                                 // it is impossible. Just the insane check
//                                 GError("Wrong boundaries. Check the code!!");
//                             } else if (prev_boundary == unispg_node->end) {
//                                 if (lclg_node->start < prev_boundary) {
// fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
// fprintf(stderr,"\n  >>>>  Graph node: |(s)----------.......(e)|\n");
// fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------....(e)|\n");
//                                     node_start_pos = prev_boundary;
//                                 } else if (lclg_node->start >= prev_boundary) {
// fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
// fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
//                                     node_start_pos = lclg_node->start;
//                                 }
//                             } else if (prev_boundary > unispg_node->end) {
//                                 // This situation has happend at least one. Insert the new unispg.
//                                 node_start_pos = lclg_node->start;
//                             }
//                             GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                             GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                             GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                             // Add new node
//                             fprintf(stderr, "^^^^ node_start_pos: %d,  unispg_node->end: %d \n", node_start_pos, lclg_node->end);
//                             node = new CGraphnodeUnispg(sample_num, node_start_pos, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                             new_no2gnode_unispg[unispg_idx].Add(node);
//                             new_unispg_nodeid += 1;
//                             // !!!!!! Do not update the prev_boundary!!!
//                             // prev_boundary = lclg_node->end;
//                             // Move to the next node of the lclg
//                             lclg_idx += 1;
//                         } else if (!lclg_is_end && !unispg_is_end) {
//                             // This is an insane check. If any of unispg or lclg is not end, at least 1 must move. 
//                                 GError("Wrong boundaries. Check the code!!");
//                         } else if (lclg_is_end && unispg_is_end) {
//                             if (unispg_node->start < lclg_node->start) {
//                                 if (unispg_node->end < lclg_node->start) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
//                                     AddBoundary(boundaries, lclg_node->start, boundaries_types, LCLG_S);
//                                     AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
                                    
//                                     GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                     GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                     GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                     // Add new node
//                                     node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                     new_no2gnode_unispg[unispg_idx].Add(node);
//                                     new_unispg_nodeid += 1;
//                                     prev_boundary = lclg_node->end;

//                                 } else if (unispg_node->end == lclg_node->start) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
//                                     AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);

//                                     GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                     GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                     GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                     // Add new node
//                                     node = new CGraphnodeUnispg(sample_num, lclg_node->start, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                     new_no2gnode_unispg[unispg_idx].Add(node);
//                                     new_unispg_nodeid += 1;
//                                     prev_boundary = lclg_node->end;

//                                 } else if (lclg_node->start < unispg_node->end) {
//                                     if (unispg_node->end < lclg_node->end) {
//                                         fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
//                                         AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                                         GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                         GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                         GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                         // Add new node
//                                         node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                         new_no2gnode_unispg[unispg_idx].Add(node);
//                                         new_unispg_nodeid += 1;
//                                         prev_boundary = lclg_node->end;

//                                     } else if (unispg_node->end == lclg_node->end) {
//                                         fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----------(e)|\n");
//                                         // Do not need to create a new node
//                                     } else if (lclg_node->end < unispg_node->end) {
//                                         fprintf(stderr,"\n  >>>>  Graph node: -----|(s)----------(e)|-----\n");
//                                         AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                         // Add new node
//                                         node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                                         new_no2gnode_unispg[unispg_idx].Add(node);
//                                         new_unispg_nodeid += 1;
//                                         prev_boundary = unispg_node->end;   
//                                     }
//                                 }
//                             } else if (unispg_node->start == lclg_node->start) {
//                                 if (unispg_node->end < lclg_node->end) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: |(s)----------.......(e)|\n");
//                                     AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                                     GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                     GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                     GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                     // Add new node
//                                     node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                     new_no2gnode_unispg[unispg_idx].Add(node);
//                                     new_unispg_nodeid += 1;
//                                     prev_boundary = lclg_node->end;

//                                 } else if (unispg_node->end == lclg_node->end) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: |(s)------------(e)|\n");
//                                 } else if (lclg_node->end < unispg_node->end) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: |(s)---------(e)|---\n");
//                                     AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                     // Add new node
//                                     node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                                     new_no2gnode_unispg[unispg_idx].Add(node);
//                                     new_unispg_nodeid += 1;
//                                     prev_boundary = unispg_node->end;   
//                                 }
//                             } else if (lclg_node->start < unispg_node->start) {
//                                 if (unispg_node->start < lclg_node->end) {
//                                     if (unispg_node->end < lclg_node->end) {
//                                         fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------....(e)|\n");
//                                         AddBoundary(boundaries, lclg_node->end, boundaries_types, LCLG_E);
//                                         GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
//                                         GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
//                                         GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
//                                         // Add new node
//                                         node = new CGraphnodeUnispg(sample_num, unispg_node->end, lclg_node->end, new_unispg_nodeid, is_passed_s, cov_s, capacity_s, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
//                                         new_no2gnode_unispg[unispg_idx].Add(node);
//                                         new_unispg_nodeid += 1;
//                                         prev_boundary = lclg_node->end;

//                                     } else if (unispg_node->end == lclg_node->end) {
//                                         fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|\n");;
//                                     } else if (lclg_node->end < unispg_node->end) {
//                                         fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
//                                         AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                         // Add new node
//                                         node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                                         new_no2gnode_unispg[unispg_idx].Add(node);
//                                         new_unispg_nodeid += 1;
//                                         prev_boundary = unispg_node->end;   
//                                     }
//                                 } else if (lclg_node->end == unispg_node->start) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|----------\n");
//                                     AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                     // Add new node
//                                     node = new CGraphnodeUnispg(sample_num, lclg_node->end, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                                     new_no2gnode_unispg[unispg_idx].Add(node);
//                                     new_unispg_nodeid += 1;
//                                     prev_boundary = unispg_node->end;  
//                                 } else if (lclg_node->end < unispg_node->start) {
//                                     fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n \n");
//                                     AddBoundary(boundaries, unispg_node->start, boundaries_types, UNISPG_S);
//                                     AddBoundary(boundaries, unispg_node->end, boundaries_types, UNISPG_E);
//                                     // Add new node
//                                     node = new CGraphnodeUnispg(sample_num, unispg_node->start, unispg_node->end, new_unispg_nodeid, unispg_node->is_passed_s, unispg_node->cov_s, unispg_node->capacity_s, false, 0, 0, 0);
//                                     new_no2gnode_unispg[unispg_idx].Add(node);
//                                     new_unispg_nodeid += 1;
//                                     prev_boundary = unispg_node->end;  
//                                 }
//                             }
//                             // if (lclg_node->start < unispg_node->start) {
//                             //     if (lclg_node->end < unispg_node->start) {
//                             //         fprintf(stderr,"\n  >>>>  Graph node: |(s)..........(e)|  ----------\n");
//                             //     } else if (lclg_node->end >= unispg_node->start) {
//                             //         fprintf(stderr,"\n  >>>>  Graph node: |(s)....----------(e)|---\n");
//                             //     }
//                             // } else if (lclg_node->start == unispg_node->start) {
//                             // } else if (lclg_node->start > unispg_node->start) {
//                             //     if (lclg_node->start >unispg_node->end) {
//                             //         fprintf(stderr,"\n  >>>>  Graph node: ----------  |(s).................(e)|\n");
//                             //     } else if (lclg_node->start == unispg_node->end) {
//                             //         fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n");
//                             //     } else if (lclg_node->start < unispg_node->end) {
//                             //         fprintf(stderr,"\n  >>>>  Graph node: ------|(s)----.............(e)|\n");
//                             //     }
//                             // }

//                             break;
//                         }
//                     }
//                 }
                
//                 GVec<bool>* is_passed_s_sink = new GVec<bool>(sample_num-1, false);
//                 GVec<float>* cov_s_sink = new GVec<float>(sample_num-1, 0.0f);
//                 GVec<float>* capacity_s_sink = new GVec<float>(sample_num-1, 0.0f);
//                 CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, new_unispg_nodeid, is_passed_s_sink, cov_s_sink, capacity_s_sink, false, 0, 0, 0);

//                 new_no2gnode_unispg[unispg_idx].Add(sink);
//                 // CGraphnodeUnispgType = UNI_LCL_NODE;

//                 no2gnode_unispg[s][unispg_idx].Clear();
//                 no2gnode_unispg[s][unispg_idx] = new_no2gnode_unispg[unispg_idx];

//                 fprintf(stderr, "Boundaries: ");
//                 for (int i = 0; i < boundaries.Count(); i++) {
//                     fprintf(stderr, "%d, ", boundaries.Get(i));
//                 }
            }



            if (more_unispg) {
                fprintf(stderr, "&&&& Add more more_unispg\n");
                current_gidx[s] += 1;
                unispg_merge_idx += 1;
            }
            if (more_lclg) {
                fprintf(stderr, "&&&& Add more more_lclg\n");

                // no2gnode = no2gnode_base + lclg_idx;

                // fprintf(stderr, ">> no2gnode: %p\n", no2gnode);
                // fprintf(stderr, ">> count: %d \n", no2gnode->Count());
                // uint lclg_start = no2gnode->Get(1)->start;
                // uint lclg_end = no2gnode->Get(no2gnode->Count()-2)->end;
                // fprintf(stderr, ">> boundary: %d - %d \n", lclg_start, lclg_end);

                lclg_idx += 1;
                lclg_merge_idx += 1;


                if (lclg_idx == lclg_limit-1) {
                    break;
                }

                // while (lclg_idx < lclg_limit-1) {
                //     lclg_idx += 1;
                //     lclg_merge_idx += 1;
                //     // Make sure the next graph is not empty
                //     int node_num = (no2gnode_base + lclg_idx)->Count();
                //     if (node_num == 0) {
                //         continue;
                //     } else {
                //         break;
                //     }
                // }


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

void UnispgGp::MoveUnispg(bool& unispg_is_end, bool& unispg_move) {
    if (!unispg_is_end) {
        unispg_move = true;
    }
}

void UnispgGp::MoveLclg(bool& lclg_is_end, bool& lclg_move) {
    if (!lclg_is_end) {
        lclg_move = true;
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