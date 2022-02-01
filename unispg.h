#ifndef __UNISPG_H__
#define __UNISPG_H__
#include "rlink.h"
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "time.h"
#include "tablemaker.h"
#include "GHashMap.hh"
#include <fstream>
#include <regex>
#include <string>
#include <typeinfo>

#include <vector>
using namespace std;


// struct Unispg {
//    protected:
//         int refstart;
//         int refend;
//         int s;
//         // int g_idx;
//         int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
//         int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
//         //  Source and sink are also included.!!
//         GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
//    public:
//         // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
//         // b: all bundles on all strands: 0,1,2
//    	    Unispg(int refstart_i=0, int refend_i=0, int s_i=0):refstart(refstart_i),refend(refend_i),s(s_i) { 
//             no2gnode = new GPVec<CGraphnode>;
//             no2gnode->setCapacity(8192);
//             // fprintf(stderr, "no2gnode print test %d: ", no2gnode);
//             CGraphnode* source=new CGraphnode(0,0,0);
//             no2gnode -> Add(source);
//             // fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(0)->nodeid);
//             graphno = 1;
//             edgeno = 0;
//         }
//         ~Unispg() {
//             // no2gnode->Clear();
//         }
//         void Clear() {
//             no2gnode -> Clear();
//             no2gnode->setCapacity(8192);
//         }
//         void AddNode(int refstart_i, int refend_i, int start, int end, int node_id, float cov) {
//             if (refstart_i == refstart && refend_i == refend) {
//             CGraphnode* new_node=new CGraphnode(start, end, node_id, cov); // start,end,nodeno
//             no2gnode -> Add(new_node);
//             // fprintf(stderr, "Node size: %d\n ", no2gnode->Count());
//             graphno += 1;
//             }
//         }
//         void AddSink(int refstart_i, int refend_i) {
//             if (refstart_i == refstart && refend_i == refend) {
//             CGraphnode* sink=new CGraphnode(0,0,graphno);
//             no2gnode -> Add(sink);
//             // fprintf(stderr, "sink Node size: %d\n ", no2gnode->Count());
//             // fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(graphno)->nodeid);
//             graphno += 1;
//             }
//         }
//         void AddEdge(int refstart_i, int refend_i, int head, int tail) {
//             if (refstart_i == refstart && refend_i == refend) {
//             // no2gnode[head];
//             // fprintf(stderr, "no2gnode[0]->nodeid %d: ", no2gnode->Get(head));
//             no2gnode->Get(head)->child.Add(tail);
//             no2gnode->Get(tail)->parent.Add(head);
//             // GBitVec childpat;
//             // GBitVec parentpat;
//             edgeno+=1;
//             }
//         }
//         void PrintGraph() {
//             { //DEBUG ONLY
//             fprintf(stderr,"\tafter traverse:\n");
//             for(int i=1;i<graphno-1;i++) {
//                 fprintf(stderr,"\tNode %d with parents:",i);
//                 for(int p=0;p<no2gnode->Get(i)->parent.Count();p++) fprintf(stderr," %d",no2gnode->Get(i)->parent[p]);
//                 fprintf(stderr," and children:");
//                 for(int c=0;c<no2gnode->Get(i)->child.Count();c++) fprintf(stderr," %d",no2gnode->Get(i)->child[c]);           
//                 fprintf(stderr," %d(%d-%d)",i,no2gnode[0][i]->start,no2gnode[0][i]->end);
//                 fprintf(stderr,"\n");
//             }
//             }
//         }
//         int get_refstart() {
//             return refstart;
//         }
//         int get_refend() {
//             return refend;
//         }
//         int get_s() {
//             return s;
//         }
//         int get_graphno() {
//             return graphno;
//         }
//         int get_edgeno() {
//             return edgeno;
//         }
//         GPVec<CGraphnode>* get_no2gnode() {
//             return no2gnode;
//         }
//     // after the graph is created, I need to do 'predict transcripts for unstranded bundles here'
//     /*****************************
//      ** 'CPrediction': constructor
//     **		predict transcripts for unstranded bundles here
//     *****************************/
//     // And then, I need to parse the graph here
//     /*****************************
//      ** 5. parse graph
//     **    'process_refguides' & 'process_transfrags' & 'find_transcripts' & 'free_treepat'
//     *****************************/
// };

struct CGraphnodeUnispg:public GSeg {
    int sample_num = 0;
	int nodeid;
	// float cov;
	// float capacity; // sum of all transcripts abundances exiting and through node
	// float rate; // conversion rate between in and out transfrags of node
	//float frag; // number of fragments included in node

// Samples vector
    // GVec<GStr> samples;

// Samples having this node
    GVec<bool> is_passed_s;
// Node coverage for each sample
    GVec<float> cov_s;
// Node capacity for each sample
    GVec<float> capacity_s;

    // GVec<>;
	GVec<int> child;
	GVec<int> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	CGraphnodeUnispg(int s=0,int e=0,unsigned int id=MAX_NODE, GVec<bool> is_passed_s_i=NULL, bool is_passed=false, GVec<float> cov_s_i=NULL, float cov=0, GVec<float> capacity_s_i=NULL, float capacity=0,float r=0):GSeg(s,e),nodeid(id),is_passed_s(),cov_s(),capacity_s(),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){
        is_passed_s = is_passed_s_i.Add(is_passed);
        cov_s = cov_s_i.Add(cov);
        capacity_s = capacity_s_i.Add(capacity);
    }

    void setup_parent() {

    }

    void setup_child() {

    }



	// CGraphnode(CGraphnode* node) {
	// 	// fprintf(stderr, "Copying node start!! //\n");
	// 	this->start = node->start;
	// 	this->end = node->end;
	// 	fprintf(stderr, "Start - end: %u - %u !!\n", this->start, this->end);
	// 	nodeid = node->nodeid;
	// 	cov = node->cov;
	// 	capacity = node->capacity;
	// 	rate = node->rate;
	// 	GVec<int>* child_cp = new GVec<int>(node->child);
	// 	child = *child_cp;
	// 	// for (int i = 0; i < node->child.Count(); i++) {
	// 	// 	child.Add(node->child.Get(i));
	// 	// }
	// 	GVec<int>* parent_cp = new GVec<int>(node->parent);
	// 	parent = *parent_cp;
	// 	// for (int i = 0; i < node->parent.Count(); i++) {
	// 	// 	parent.Add(node->parent.Get(i));
	// 	// }

		
	// 	// childpat = node->childpat;
	// 	// parentpat = node->parentpat;
	// 	GVec<int>* trf_cp = new GVec<int>(node->trf);
	// 	trf = *trf_cp;
	// 	// for (int i = 0; i < node->trf.Count(); i++) {
	// 	// 	trf.Add(node->trf.Get(i));
	// 	// }
	// 	hardstart = node->hardstart;
	// 	hardend = node->hardend;
	// 	// fprintf(stderr, "End of copying node!!\n");
	// }
};


struct UnispgGp {
    protected:
        GPVec<CGraphnodeUnispg>* no2gnode_unispg[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
        GVec<int> current_gidx;
        GVec<GStr> samples;
    public:
        // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
        // b: all bundles on all strands: 0,1,2
        UnispgGp() { 
            for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
                int s=sno/2; // adjusted strand due to ignoring neutral strand
                no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
			    // current_gidx[s] = 0;
            }
        }

        void ProcessSample(GStr sample_name) {
            samples.Add(sample_name);
            for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
                int s=sno/2; // adjusted strand due to ignoring neutral strand
			    current_gidx[s] = 0;
            }
        }

        ~UnispgGp() {
            for(int i=0;i<2;i++) {
            delete [] no2gnode_unispg[i];
            };
        }

        // void SetUnispgCapacity(int s) {
        //     fprintf(stderr, "**************************************\n");
        //     fprintf(stderr, "*********** SetUnispgCapacity ********\n");
        //     fprintf(stderr, "**************************************\n");
        //     no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
        // }

        void AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode) {   
            if (fidx == 0) {
                fprintf(stderr, "*****************************\n");
                fprintf(stderr, "*********** AddGraph ********\n");
                fprintf(stderr, "*****************************\n");
                int cgidx = current_gidx[s];
                // This is the first unispg. Simply add copy it into UnispgGp
                for (int i=0; i<no2gnode->Count(); i++) {
                    CGraphnode* node = new CGraphnode(no2gnode->Get(i));

                    GVec<bool> is_passed_s = NULL;
                    GVec<float> cov_s = NULL;
                    GVec<float> capacity_s = NULL;

                    CGraphnodeUnispg* node_unispg = new CGraphnodeUnispg(node->start, node->end, node->nodeid, is_passed_s, 1, cov_s, node->cov, capacity_s, node->capacity, node->rate);

                    // no2gnodeGp_unispg[s][cgidx].Add(node);
                    no2gnode_unispg[s][cgidx].Add(node_unispg);

                    // fprintf(stderr, "Inside!! cgidx: %d \n", cgidx);
                    // fprintf(stderr, "no2gnode->Get(i)->nodeid: %d \n", node->nodeid);
                    // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx].Last()->nodeid: %d \n", no2gnodeGp_unispg[s][cgidx].Last()->nodeid);
                    // fprintf(stderr, "no2gnodeGp_unispg[s][cgidx].Last() start - end : %d - %d \n", int(no2gnodeGp_unispg[s][cgidx].Last()->start), int(no2gnodeGp_unispg[s][cgidx].Last()->end));
                    // fprintf(stderr, "no2gnode_unispg[s][cgidx].Last()->nodeid: %d \n", no2gnode_unispg[s][cgidx].Last()->nodeid);
                    // fprintf(stderr, "no2gnode_unispg[s][cgidx].Last() start - end : %d - %d \n", int(no2gnode_unispg[s][cgidx].Last()->start), int(no2gnode_unispg[s][cgidx].Last()->end));
                }

                // Linking parent and child
                for (int i=0; i<no2gnode->Count(); i++) {
                    CGraphnode* node = new CGraphnode(no2gnode->Get(i));
                    fprintf(stderr, "CGraphnode parent: ");
                    for(int p=0; p<no2gnode->Get(i)->parent.Count(); p++) {
                        fprintf(stderr, "%d  ", no2gnode->Get(i)->parent[p]);
                        fprintf(stderr, "%d  ", no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->parent[p] ]->nodeid);

                        no2gnode_unispg[s][cgidx][i]->parent.Add(no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->parent[p] ]->nodeid);
                    }
                    fprintf(stderr, "\nCGraphnode child: ");
                    for(int c=0; c<no2gnode->Get(i)->child.Count(); c++) {
                        fprintf(stderr, "%d  ", no2gnode->Get(i)->child[c]); 
                        fprintf(stderr, "%d  ", no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->child[c] ]->nodeid); 

                        no2gnode_unispg[s][cgidx][i]->child.Add(no2gnode_unispg[s][cgidx][ no2gnode->Get(i)->child[c] ]->nodeid);
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
                int graph_overlap = false;
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
                        graph_overlap = false;
                    } else if (unispg_start < lclg_start && unispg_end >= lclg_start && unispg_end <= lclg_end) {
                        // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
                        fprintf(stderr,"\n  &&& Graph: ----------|(s).................(e)|   or   -----|(s)-----............(e)| \n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start < lclg_start && unispg_end > lclg_end) {
                        // -----|(s)------------(e)|--
                        fprintf(stderr,"\n &&& Graph: -----|(s)------------(e)|-- \n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start == lclg_start && unispg_end < lclg_end) {
                        // |(s)----------.................(e)| 
                        fprintf(stderr,"\n &&& Graph: |(s)----------............(e)|\n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start > lclg_start && unispg_end < lclg_end) {
                        // |(s)........----------........(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)........----------........(e)| \n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start > lclg_start && unispg_end == lclg_end) {
                        // |(s)............----------(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)............----------(e)| \n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start == lclg_start && unispg_end == lclg_end) {
                        // |(s)----------(e)|
                        fprintf(stderr,"\n &&& Graph: |(s)----------(e)| \n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start <= lclg_end && unispg_end > lclg_end) {
                        // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
                        fprintf(stderr,"\n &&& Graph: (s)...............------(e)|-----    or   |(s).................(e)|---------- \n");
                        next_unispg = true;
                        graph_overlap = true;
                    } else if (unispg_start > lclg_end) {
                        // The node is outside the current bundle => This node belongs to the next bundlenode
                        // |(s).................(e)|   ----------
                        fprintf(stderr,"\n &&& Graph: |(s).................(e)|   ---------- \n");
                        next_unispg = false;
                        graph_overlap = false;
                        current_gidx[s] -= 1;
                    } else {
                        fprintf(stderr,"\n &&& Unknown area!!!! \n");
                        next_unispg = false;
                        graph_overlap = false;
                    }
                    for (int i = 0; i < no2gnode->Count(); i++) {
                        // no2gnode->Get(i);
                        fprintf(stderr, "&&&& This is the local graphnode: %d\n", no2gnode->Get(i)->nodeid);
                    }


                    /*****************************
                     * Create a new graph 
                     *****************************/
                    if (graph_overlap == true) {

                        unsigned int nodeid = 0;
                        // no2gnode_unispg[s][cgidx];

                        GPVec<CGraphnodeUnispg>* new_no2gnode_unispg; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
                        new_no2gnode_unispg = new GPVec<CGraphnodeUnispg>[20000];


                        int lclg_idx = 1;
                        int unispg_idx = 1;
                        bool more_comparison = true;
                        bool lclg_is_end = false;
                        bool unispg_is_end = false;
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
            // CGraphnodeUnispg(int s=0,int e=0,unsigned int id=MAX_NODE, GVec<bool> is_passed_s_i=NULL, bool is_passed=false, GVec<float> cov_s_i=NULL, float cov=0, GVec<float> capacity_s_i=NULL, float capacity=0,float r=0):GSeg(s,e),nodeid(id),is_passed_s(),cov_s(),capacity_s(),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){

        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // To create a universal graph node, I should input a universal graphnode as well as a local graphnode.
        CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(unispg_node->start, unispg_node->end, nodeid, unispg_node->is_passed_s, false, unispg_node->cov_s, 0, unispg_node->capacity_s, 0, 0);
        CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(lclg_node->start, lclg_node->end, nodeid, NULL, true, NULL, lclg_node->cov, NULL, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);


                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                            } else if (unispg_node->start < lclg_node->start && unispg_node->end == lclg_node->start) {
                                fprintf(stderr,"\n  >>>>  Graph node: ----------|(s).................(e)|\n\n");
        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        CGraphnodeUnispg* node_1 = new CGraphnodeUnispg(unispg_node->start, unispg_node->end, nodeid, false, 0, 0, 0);
        CGraphnodeUnispg* node_2 = new CGraphnodeUnispg(lclg_node->start, lclg_node->end, nodeid, true, lclg_node->cov, lclg_node->capacity, lclg_node->rate);

                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                            } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->start && unispg_node->end < lclg_node->end) {
                                fprintf(stderr,"\n  >>>>  Graph node: -----|(s)-----............(e)| \n\n");

        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, lclg_node->start, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(lclg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_3 = new CGraphnode(unispg_node->end, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                            } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->start && unispg_node->end == lclg_node->end) {
                                fprintf(stderr,"\n  >>>>  Graph node: -----|(s)---------(e)| \n\n");

        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, lclg_node->start, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);

                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start < lclg_node->start && unispg_node->end > lclg_node->end) {
                                fprintf(stderr,"\n >>>>  Graph node: -----|(s)------------(e)|-- \n\n");


            // CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0):GSeg(s,e),
            // 		nodeid(id),cov(nodecov),capacity(cap),rate(r),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){}

        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, lclg_node->start, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_3 = new CGraphnode(lclg_node->end, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);


                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start == lclg_node->start && unispg_node->end < lclg_node->end) {
        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->end, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
                                fprintf(stderr,"\n >>>>  Graph node: |(s)----------............(e)|\n\n");
                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                            } else if (unispg_node->start == lclg_node->start && unispg_node->end == lclg_node->end) {
        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
                                fprintf(stderr,"\n >>>>  Graph node: |(s)----------(e)| \n\n");
                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start == lclg_node->start && unispg_node->end > lclg_node->end) {
                                fprintf(stderr,"\n >>>>  Graph node: |(s)----------(e)|------ \n\n");
        // CGraphnode* node_1 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->end, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start > lclg_node->start && unispg_node->end < lclg_node->end) {
                                fprintf(stderr,"\n >>>>  Graph node: |(s)........----------........(e)| \n\n");
        // CGraphnode* node_1 = new CGraphnode(lclg_node->start, unispg_node->start, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_3 = new CGraphnode(unispg_node->end, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                            } else if (unispg_node->start > lclg_node->start && unispg_node->end == lclg_node->end) {
                                fprintf(stderr,"\n >>>>  Graph node: |(s)............----------(e)| \n\n");
        // CGraphnode* node_1 = new CGraphnode(lclg_node->start, unispg_node->start, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);

                                if (!unispg_is_end) {
                                    unispg_move = true;
                                }
                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start < lclg_node->start && unispg_node->start > lclg_node->end && unispg_node->end > lclg_node->end) {
                                fprintf(stderr,"\n >>>>  Graph node: |(s)...............------(e)|----- \n\n");
        // CGraphnode* node_1 = new CGraphnode(lclg_node->start, unispg_node->start, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->start, lclg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
        // CGraphnode* node_3 = new CGraphnode(lclg_node->end, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start == lclg_node->end && unispg_node->end > lclg_node->end) {
                                fprintf(stderr,"\n >>>>  Graph node: |(s).................(e)|---------- \n\n");
        // CGraphnode* node_1 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
                                if (!lclg_is_end) {
                                    lclg_move = true;
                                }
                            } else if (unispg_node->start > lclg_node->end) {
                                // The node is outside the current bundle => This node belongs to the next bundlenode
                                fprintf(stderr,"\n >>>>  Graph node: |(s).................(e)|   ---------- \n\n");
        // CGraphnode* node_1 = new CGraphnode(lclg_node->start, lclg_node->end, 0, lclg_node->cov, lclg_node->capacity, lclg_node->rate);
        // CGraphnode* node_2 = new CGraphnode(unispg_node->start, unispg_node->end, 0, unispg_node->cov, unispg_node->capacity, unispg_node->rate);
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
                        }
                    }
                    current_gidx[s] += 1;
                }



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
                delete [] no2gnode_unispg[i];
                no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
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

        // int get_refstart () {
        //     return refstart;
        // }
        // int get_refend () {
        //     return refend;
        // }
        // GVec<int>* get_graphnoGp () {
        //     return graphnoGp;
        // }
        // GVec<int>* get_edgenoGp () {
        //     return edgenoGp;
        // }
        GPVec<CGraphnodeUnispg>** get_no2gnodeGp () {
            return no2gnode_unispg;
        }
};







// track_idx = 0;
// for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
//     int s=sno/2; // adjusted strand due to ignoring neutral strand
// //     gpSize[s] = 0;
// //     graphnoGp[s].setCapacity(8192);
// //     edgenoGp[s].setCapacity(8192);
// //     no2gnodeGp[s] = new GPVec<CGraphnode>[8192];
//     no2gnodeGp_unispg[s] = new GPVec<CGraphnode>[20000];
// }



// void unispgGp_clean();

// void AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode);

// void Clear();
// // Need to be modifed!!!
// void PrintGraphGp();


#endif