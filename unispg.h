#ifndef __UNISPG_H__
#define __UNISPG_H__
#pragma once
#include "global_params.h"
// #include "rlink.h"

// extern int sample_num;
// extern GVec<FILE*> pos_dot_vec;
// extern GVec<GStr> pos_dotfname_vec; 
// extern GVec<FILE*> neg_dot_vec;
// extern GVec<GStr> neg_dotfname_vec; 

// extern GVec<FILE*>* dots[2];
// extern GVec<GStr>* dotfname_vec[2];


// extern GVec<FILE*>* node_lclg_bed_vec[2];
// extern GVec<GStr>* nodelclgfname_vec[2]; 
// extern GVec<FILE*>* edge_lclg_bed_vec[2];
// extern GVec<GStr>* edgelclgfname_vec[2]; 

// extern GVec<FILE*>* node_novp_bed_vec[2];
// extern GVec<GStr>* nodenovpfname_vec[2]; 
// extern GVec<FILE*>* edge_novp_bed_vec[2];
// extern GVec<GStr>* edgenovpfname_vec[2]; 

// extern GVec<FILE*>* node_unispg_bed_vec[2];
// extern GVec<GStr>* nodeunispgfname_vec[2]; 
// extern GVec<FILE*>* edge_unispg_bed_vec[2];
// extern GVec<GStr>* edgeunispgfname_vec[2]; 
// extern FILE* node_unispg_unstrand_bed; 


typedef std::tuple<int, int, int> b_g_n_tuple;

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& p) const {
		return std::hash<T1>{}(p.first) ^ std::hash<T2>{}(p.second);
		// return std::hash<T1>{}(std::get<0>(p)) ^ std::hash<T2>{}(std::get<1>(p)) ^ std::hash<T3>{}(std::get<2>(p));
	}
};

struct tuple_hash {
	template <class T1, class T2, class T3>
	std::size_t operator() (const std::tuple<T1, T2, T3>& p) const {
		return std::hash<T1>{}(std::get<0>(p)) ^ std::hash<T2>{}(std::get<1>(p)) ^ std::hash<T3>{}(std::get<2>(p));
		// return std::hash<T1>{}(std::get<0>(p)) ^ std::hash<T2>{}(std::get<1>(p)) ^ std::hash<T3>{}(std::get<2>(p));
	}
};

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

struct CGraphnodeUnispg:public GSeg {
    int sample_num = 0;
	int old_graph_id = -1;
	int old_node_id = -1;

	/*****************************
	 ** APPLY_UNISPG & CREATE_UNISPG
	 *****************************/
	int nodeid;
	// Samples having this node
    GVec<bool> is_passed_s;
	// Node coverage for each sample
    GVec<float> cov_s;
	// Node capacity for each sample
    GVec<float> capacity_s;
    GVec<float> rate_s;
	GVec<int> child;
	GVec<int> parent;

	// Node coverage for each sample
    // GPVec<float> cov_unispg_s;
	float* cov_unispg_s;
	// [2] = {0};  // how many nodes are in a certain graph g, on strand s: graphno[s][g]


	/*****************************
	 ** To-do 
	 *****************************/
	GBitVec childpat;
	GBitVec parentpat;
	GVec<int> trf; // transfrags that pass the node
	bool hardstart:1; // verified/strong start
	bool hardend:1;	// verified/strong end
	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
	// CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0, int old_graph_id_i=0, int old_node_id_i=0, unsigned int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, GVec<float>* cov_s_i=NULL, GVec<float>* capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0):GSeg(s,e),sample_num(sample_num_i), old_graph_id(old_graph_id_i), old_node_id(old_node_id_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){
	CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0, int id=MAX_NODE, GVec<bool> is_passed_s_i=NULL, GVec<float> cov_s_i=NULL, GVec<float> capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0, bool set_g_n_idx=false, int g_idx=-1, int n_idx=-1):GSeg(s,e),sample_num(sample_num_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i), child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){

		// fprintf(stderr, "		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		// fprintf(stderr, "		^^^ Creating graphnode id: %d (%u - %u) \n", id, s, e);
		// fprintf(stderr, "		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		// fprintf(stderr, "		^^ is_passed: %d\n", is_passed);
		// fprintf(stderr, "		^^ cov      : %f\n", cov);
		// fprintf(stderr, "		^^ capacity : %f\n", capacity);

		// if (is_passed_s_i == NULL) {
        //     GVec<bool> is_passed_source;
		// }
		// if (cov_s_i == NULL) {
        //     GVec<float> cov_source;
		// }
		// if (capacity_s_i == NULL) {
        //     GVec<float> capacity_source;
		// }
		// fprintf(stderr, "sample_num: %d\n", sample_num);
		is_passed_s.cAdd(is_passed);
		// fprintf(stderr, "is_passed_s[0]: %d\n", is_passed_s->Get(0));

		cov_s.cAdd(cov);
		// fprintf(stderr, "cov_s[0]: %F\n", cov_s->Get(0));

		capacity_s.cAdd(capacity);
		// fprintf(stderr, "capacity_s[0]: %f\n", capacity_s->Get(0));

		// To-do: add it as a variable that need to be passed into the constructor.
        // GVec<float> rate_s;
		//  = new GVec<float>(sample_num-1, 0.0f);

		// fprintf(stderr, "rate: %f\n", r);
		rate_s.cAdd(r);

		// fprintf(stderr, "sample_num: %d\n", sample_num);
		cov_unispg_s = new float[sample_num];

		if (set_g_n_idx) {
			old_graph_id = g_idx;
			old_node_id = n_idx;
		}
    }

	void set_nodeid(int new_node_id) {
		nodeid = new_node_id;
	}

	// void add_cov_unispg_s(float cov) {
	// 	cov_unispg_s->cAdd(cov);
	// }

	int calOverlapLen(uint rstart, uint rend) {
		// fprintf(stderr, ">>>> Inside calOverlapLen:!!!\n");
		// fprintf(stderr, ">> (%u - %u); (%u - %u)!!!\n", rstart, rend, start, end);

		if (rstart>rend) { Gswap(rstart,rend); }
		if (start<rstart) {
			if (rstart>end) return 0;
			return (rend>end) ? end-rstart+1 : rend-rstart+1;
		}
		else { //rstart<=start
			if (start>rend) return 0;
			return (rend<end)? rend-start+1 : end-start+1;
		}
	}
};

struct UnispgGp {
	public:
		int refstart = 0; // the start of the first node.
		int refend = 0; // the end of the last node.
		GStr refseq;
		GPVec<CGraphnodeUnispg>* no2gnode_unispg[2]; // for each graph g, on a strand s, no2gnode_unispg[s][g][i] gives the node i

		int graph_num[2] = {0};  // how many graph does a unispgs has.
		GVec<int> node_nums[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
		GVec<int> edge_nums[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]

		GVec<GStr> samples;
		CGraphnodeUnispg* source_gp[2];
		CGraphnodeUnispg* sink_gp[2];

		UnispgGp() {
			for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
				int s=sno/2; // adjusted strand due to ignoring neutral strand
				no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
			}
		}

		~UnispgGp() {
			for(int i=0;i<2;i++) {
				delete [] no2gnode_unispg[i];
			};
		}

		int get_refstart() {
			return refstart;
		}

		int get_refend() {
			return refend;
		}

		void Clear_no2gnode_unispg() {
			for(int i=0;i<2;i++) {
				delete [] no2gnode_unispg[i];
				no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
		}

		void ProcessSamples(GVec<GStr> sample_nams);
		void ProcessSample(GStr sample_name);
		void PrintGraphGp();
		GPVec<CGraphnodeUnispg>** get_no2gnodeGp ();
};
#endif

