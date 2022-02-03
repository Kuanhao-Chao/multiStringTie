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
    GVec<bool>* is_passed_s;
// Node coverage for each sample
    GVec<float>* cov_s;
// Node capacity for each sample
    GVec<float>* capacity_s;

    // GVec<>;
	GVec<CGraphnodeUnispg*> child;
	GVec<CGraphnodeUnispg*> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0,unsigned int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, GVec<float>* cov_s_i=NULL, GVec<float>* capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0):GSeg(s,e),sample_num(sample_num_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){
// for (int i = 0; i < is_passed_s_i->Count(); i++) {
//     fprintf(stderr, "is_passed_s_i: %d\n", is_passed_s_i->Get(i));
// }
// for (int i = 0; i < cov_s_i->Count(); i++) {
//     fprintf(stderr, "cov_s_i: %f\n", cov_s_i->Get(i));
// }
// for (int i = 0; i < capacity_s_i->Count(); i++) {
//     fprintf(stderr, "capacity_s_i: %f\n", capacity_s_i->Get(i));
// }

// for (int i = 0; i < is_passed_s->Count(); i++) {
//     fprintf(stderr, "is_passed_s: %d\n", is_passed_s->Get(i));
// }
// for (int i = 0; i < cov_s->Count(); i++) {
//     fprintf(stderr, "cov_s: %f\n", cov_s->Get(i));
// }
// for (int i = 0; i < capacity_s->Count(); i++) {
//     fprintf(stderr, "capacity_s: %f\n", capacity_s->Get(i));
// }
		is_passed_s->cAdd(is_passed);
		cov_s->cAdd(cov);
		capacity_s->cAdd(capacity);			
// for (int i = 0; i < is_passed_s->Count(); i++) {
//     fprintf(stderr, "ADD is_passed_s: %d\n", is_passed_s->Get(i));
// }
// for (int i = 0; i < cov_s->Count(); i++) {
//     fprintf(stderr, "ADD cov_s: %f\n", cov_s->Get(i));
// }
// for (int i = 0; i < capacity_s->Count(); i++) {
//     fprintf(stderr, "ADD capacity_s: %f\n", capacity_s->Get(i));
// }
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

        ~UnispgGp() {
            for(int i=0;i<2;i++) {
            delete [] no2gnode_unispg[i];
            };
        }
        
        void ProcessSample(GStr sample_name);

        void AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode);

        void Clear() {
            // fprintf(stderr, "**** Start Clearing !!!! \n ");
            for(int i=0;i<2;i++) {
                delete [] no2gnode_unispg[i];
                no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
            };
        }

        void PrintGraphGp();
        GPVec<CGraphnodeUnispg>** get_no2gnodeGp ();
};


#endif