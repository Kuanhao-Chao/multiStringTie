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

extern GVec<FILE*> node_cov_pos_bed_vec;
extern GVec<FILE*> node_cov_neg_bed_vec;
extern GVec<FILE*> edge_cov_pos_bed_vec;
extern GVec<FILE*> edge_cov_neg_bed_vec;

extern GVec<FILE*> node_cov_pos_bed_unispg_vec;
extern GVec<FILE*> node_cov_neg_bed_unispg_vec;
extern GVec<FILE*> edge_cov_pos_bed_unispg_vec;
extern GVec<FILE*> edge_cov_neg_bed_unispg_vec;

extern GVec<FILE*> node_cov_pos_novp_bed_vec;
extern GVec<FILE*> node_cov_neg_novp_bed_vec;

enum LCLG_ITR_STATUS {
	OUT_OF_RANGE=0,
	LASTG_COUNT_0,
	LASTG_COUNT_N_0,
	N_LASTG_COUNT_0,
	N_LASTG_COUNT_N_0,
	LCLG_ITR_INIT
};

enum CGraphBoundaryType {
	EMPTY_TYPE=0,
	UNISPG_S,
	UNISPG_E,
	LCLG_S,
	LCLG_E,
	UNISPG_S_LCLG_S,
	UNISPG_S_LCLG_E,
	UNISPG_E_LCLG_S,
	UNISPG_E_LCLG_E
	// UNISPG_S_EMPTY=0,
	// UNISPG_S_LCLG,
	// UNISPG_E_EMPTY,
	// UNISPG_E_LCLG,
	// LCLG_S_EMPTY,
	// LCLG_S_UNISPG,
	// LCLG_E_EMPTY,
	// LCLG_E_UNISPG,
	// UNISPG_S_LCLG_S,
	// UNISPG_S_LCLG_E,
	// UNISPG_E_LCLG_S,
	// UNISPG_E_LCLG_E
};

static const char *enum_str[] = { 	
	"EMPTY_TYPE",
	"UNISPG_S",
	"UNISPG_E",
	"LCLG_S",
	"LCLG_E",
	"UNISPG_S_LCLG_S",
	"UNISPG_S_LCLG_E",
	"UNISPG_E_LCLG_S",
	"UNISPG_E_LCLG_E" 
};


enum CGraphnodeUnispgType {
	UNI_LCL_NODE=0,
	UNI_NODE,
	LCL_NODE,
	EMPTY_NODE
};

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
	int old_graph_id;
	int old_node_id;
	int nodeid;
	// Samples having this node
    GVec<bool>* is_passed_s;
	// Node coverage for each sample
    GVec<float>* cov_s;
	// Node capacity for each sample
    GVec<float>* capacity_s;

	GVec<CGraphnodeUnispg*> child;
	GVec<CGraphnodeUnispg*> parent;
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	// CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0, int old_graph_id_i=0, int old_node_id_i=0, unsigned int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, GVec<float>* cov_s_i=NULL, GVec<float>* capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0):GSeg(s,e),sample_num(sample_num_i), old_graph_id(old_graph_id_i), old_node_id(old_node_id_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){
	CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0, int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, GVec<float>* cov_s_i=NULL, GVec<float>* capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0):GSeg(s,e),sample_num(sample_num_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){

		fprintf(stderr, "		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		fprintf(stderr, "		^^^ Creating graphnode id: %d (%u - %u) \n", id, s, e);
		fprintf(stderr, "		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		fprintf(stderr, "		^^ is_passed: %d\n", is_passed);
		fprintf(stderr, "		^^ cov      : %f\n", cov);
		fprintf(stderr, "		^^ capacity : %f\n", capacity);
		is_passed_s->cAdd(is_passed);
		cov_s->cAdd(cov);
		capacity_s->cAdd(capacity);			
    }

    void setup_parent() {
    }

    void setup_child() {
    }
};


struct UnispgGp {
    public:
        GPVec<CGraphnodeUnispg>* no2gnode_unispg[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
        GVec<int> current_gidx; // graph id
        GVec<int> last_nidx; // node id
		GVec<uint> prev_bdy;
        GVec<GStr> samples;
		GVec<bool> has_unispg_tail;
		GVec<uint> new_unispg_nodeid;
		GPVec<CGraphnodeUnispg>* lclg_nonoverlap[2];

		GPVec<CGraphnodeUnispg>* new_no2gnode_unispg[2]; // for each graph g, on a strand s, no2gnode[g][i] gives the node i
        // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
        // b: all bundles on all strands: 0,1,2
        UnispgGp() { 
            for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
                int s=sno/2; // adjusted strand due to ignoring neutral strand
                no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
                new_no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
            }
        }
        ~UnispgGp() {
            for(int i=0;i<2;i++) {
            	delete [] no2gnode_unispg[i];
            	delete [] lclg_nonoverlap[i];
            	delete [] new_no2gnode_unispg[i];
            };
        }
        void ProcessSample(GStr sample_name);
		void WriteLCLG(int fidx, int s, GPVec<CGraphnode>* no2gnode, int g);
		void WriteUNISPG(int fidx, int s, int unispg_start_idx, int unispg_end_idx);
		void MergeLCLG(int s, int sample_num, GPVec<CGraphnode>* no2gnode, int lclg_limit, int boudleGP_start_idx, int boudleGP_end_idx, int& new_nonolp_lclg_idx, bool write_unispg);
		void FirstUnispgAlgo(int fidx, int s, int sample_num, GPVec<CGraphnode>* no2gnode, int lclg_limit, int& new_nonolp_lclg_idx, bool write_unispg);
		bool RecruitMRGGP(int s, int& lclg_idx, int& new_nonolp_lclg_idx, LCLG_ITR_STATUS& lclg_itr_status, bool& more_lclg, bool& try_more_unispg, int& process_ovp_graphs, int& lclg_idx_start, int& lclg_idx_end, int& unispg_idx_start, int& unispg_idx_end, int& unispg_node_idx);

		void CmpLclgNodeUnispgNode(int fidx, int s, int sample_num, CGraphnodeUnispg*& node, bool& lclg_node_move, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& unispg_node_move, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs);

		void SecondUnispgAlgo(int fidx, int s, int sample_num, CGraphnodeUnispg*& node, bool& lclg_node_move, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& lclg_next, bool& unispg_node_move, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs, bool& unispg_next);

		bool PairLclgUnispgInMRGGP(int fidx, int s, int sample_num, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& lclg_next, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs, bool& unispg_next);

		void ThirdUnispgAlgo(int fidx, int s, int sample_num, bool& lclg_reached_end, CGraphnodeUnispg*& node, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& lclg_next, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs, bool& unispg_next);

		void AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode_base, int lclg_limit);


		void WriteNonOVP(int fidx, int s, int unispg_start_idx, int unispg_end_idx);

		void MoveUnispgNode(bool& unispg_is_end, bool& unispg_move);
		void MoveLclgNode(bool& lclg_is_end, bool& lclg_move);
		void WriteGraphGp();

        void Clear() {
            // fprintf(stderr, "**** Start Clearing !!!! \n ");
            for(int i=0;i<2;i++) {
                delete [] no2gnode_unispg[i];
                no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
                delete [] new_no2gnode_unispg[i];
                new_no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
            };
        }

        void PrintGraphGp();
        GPVec<CGraphnodeUnispg>** get_no2gnodeGp ();
};


#endif