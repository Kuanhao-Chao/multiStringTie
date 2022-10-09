#ifndef __UNISPG_C_H__
#define __UNISPG_C_H__
#pragma once
#include "global_params_C.h"
#include "../unispg.h"

#include <unordered_map>
#include <iostream>

// extern UnispgGp_CREATE* unispg_gp;

struct UnispgGp_CREATE:public UnispgGp {
	public:

		/*****************************
		 ** APPLY_UNISPG & CREATE_UNISPG
		 *****************************/
		// GPVec<CTransfrag> *transfrag_unispg[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
		// CTreePat **tr2no_unispg[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
		// GIntHash<int> *gpos_unispg[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
		// GVec<int> lastgpos_unispg[2];
		GPVec<CGraphnodeUnispg>* lclg_nonoverlap[2]; // for each graph g, on a strand s, lclg_nonoverlap[s][g][i] gives the node i
		GPVec<CTransfrag>* lclg_nonoverlap_transfrag[2]; // for each transfrag t on a strand s, in a graph g, lclg_nonoverlap_transfrag[s][g][t] gives it's abundance and it's pattern

		GPVec<CGraphnodeUnispg>* new_no2gnode_unispg[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
		GVec<int> new_gidx; // graph id

		GVec<int> current_gidx; // graph id
		GVec<int> last_nidx; // node id
		GVec<uint> prev_bdy;
		GVec<int> lclg_bundle_num;
		GVec<bool> has_unispg_tail;
		GVec<int> new_unispg_nodeid;

		std::unordered_map<std::tuple<int, int, int>, GVec<int>, tuple_hash> lclg_nidx_2_new_nidx_ls_pos;
		std::unordered_map<std::tuple<int, int, int>, GVec<int>, tuple_hash> unispg_nidx_2_new_nidx_ls_pos;
		std::unordered_map<std::tuple<int, int, int>, GVec<int>, tuple_hash> lclg_nidx_2_new_nidx_ls_neg;
		std::unordered_map<std::tuple<int, int, int>, GVec<int>, tuple_hash> unispg_nidx_2_new_nidx_ls_neg;

		UnispgGp_CREATE() {
			for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
				int s=sno/2; // adjusted strand due to ignoring neutral strand
				no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
				new_no2gnode_unispg[s] = new GPVec<CGraphnodeUnispg>[20000];
				lclg_nonoverlap[s] = new GPVec<CGraphnodeUnispg>[20000];
				lclg_nonoverlap_transfrag[s] = new GPVec<CTransfrag>[20000];

				// transfrag_unispg[s] = new GPVec<CTransfrag>[20000];
				// source_gp[s] = new CGraphnodeUnispg[1];
				// sink_gp[s] = new CGraphnodeUnispg[1];
			}
		}

		~UnispgGp_CREATE() {
			for(int i=0;i<2;i++) {
				// graphno_unispg[s] = 0;
				// edgeno_unispg[s] = 0;
				delete [] no2gnode_unispg[i];
				delete [] new_no2gnode_unispg[i];
				delete [] lclg_nonoverlap[i];
				delete [] lclg_nonoverlap_transfrag[i];
				// delete [] source_gp[i];
				// delete [] sink_gp[i];
				// delete [] transfrag_unispg[i];
			};
		}

		int get_refstart() {
			return refstart;
		}

		int get_refend() {
			return refend;
		}

		void ProcessSample(GStr sample_name);
		
		void WriteLCLG(int fidx, int s, GPVec<CGraphnode>* no2gnode, int g);
		void WriteUNISPG(int fidx, int s, int unispg_start_idx, int unispg_end_idx);
		void MergeLCLG(int s, int sample_num, GPVec<CGraphnode>* no2gnode, int lclg_limit, int boudleGP_start_idx, int boudleGP_end_idx, int& new_nonolp_lclg_idx, bool write_unispg);
		void FirstUnispgAlgo(int fidx, int s, int sample_num, GPVec<CGraphnode>* no2gnode, int lclg_limit, int& new_nonolp_lclg_idx, bool write_unispg);

		bool RecruitMRGGP(int s, int& lclg_idx, int& new_nonolp_lclg_idx, LCLG_ITR_STATUS& lclg_itr_status, bool& more_lclg, bool& try_more_unispg, int& process_ovp_graphs, int& lclg_idx_start, int& lclg_idx_end, int& unispg_idx_start, int& unispg_idx_end, int& unispg_node_idx);

		int CreateThirdAlgHashEnt(int s, bool lclg_or_unispg, int g_idx, int n_idx, int new_unispg_n_idx);

		void AddThirdAlgParentEdgeSub(int s, bool lclg_or_unispg, CGraphnodeUnispg*& node, int g_i, CGraphnodeUnispg* g_node, int pre_parent);
		void AddThirdAlgParentEdge(int s, CGraphnodeUnispg*& node, int lclg_i, CGraphnodeUnispg* lclg_node, int unispg_i, CGraphnodeUnispg* unispg_node, int pre_parent_lclg, int pre_parent_unispg);
		void CmpLclgNodeUnispgNode(int fidx, int s, int sample_num, CGraphnodeUnispg*& node, bool& lclg_node_move, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& unispg_node_move, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs);

		void SecondUnispgAlgo(int fidx, int s, int sample_num, CGraphnodeUnispg*& node, bool& lclg_node_move, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& lclg_next, bool& unispg_node_move, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs, bool& unispg_next);

		bool PairLclgUnispgInMRGGP(int fidx, int s, int sample_num, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& lclg_next, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs, bool& unispg_next);

		void ThirdUnispgAlgo(int fidx, int s, int sample_num, bool& lclg_reached_end, CGraphnodeUnispg*& node, int& lclg_i, int& lclg_idx_start, int& lclg_idx_end, CGraphnodeUnispg*& lclg_node, int& lclg_node_idx, bool& lclg_is_lastnode, uint& lclg_start_pcs, uint& lclg_end_pcs, bool& lclg_next, int& unispg_i, int& unispg_idx_start, int& unispg_idx_end, CGraphnodeUnispg*& unispg_node, int& unispg_node_idx, bool& unispg_is_lastnode, uint& unispg_start_pcs, uint& unispg_end_pcs, bool& unispg_next);

		void AddGraph(int fidx, int s, GPVec<CGraphnode>* no2gnode, int lclg_limit);

		// void construct_transfrag_unispg(int fidx, int s);

		void WriteNonOVP(int fidx, int s, int unispg_start_idx, int unispg_end_idx);

		void MoveUnispgNode(bool& unispg_is_end, bool& unispg_move);
		void MoveLclgNode(bool& lclg_is_end, bool& lclg_move);
		void WriteUNISPG_DOT(int fidx, int s, int unispg_start_idx, int unispg_end_idx);

		void Clear() {
		// fprintf(stderr, "**** Start Clearing !!!! \n ");
			for(int i=0;i<2;i++) {
				delete [] no2gnode_unispg[i];
				no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];

				delete [] lclg_nonoverlap[i];
				lclg_nonoverlap[i] = new GPVec<CGraphnodeUnispg>[20000];

				delete [] lclg_nonoverlap_transfrag[i];
				lclg_nonoverlap_transfrag[i] = new GPVec<CTransfrag>[20000];

				delete [] new_no2gnode_unispg[i];
				new_no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
		}

		void Clear_lclg_nonoverlap() {
		// fprintf(stderr, "**** Start Clearing !!!! \n ");
			for(int i=0;i<2;i++) {
				delete [] lclg_nonoverlap[i];
				lclg_nonoverlap[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
		}

		void Clear_lclg_nonoverlap_transfrag() {
		// fprintf(stderr, "**** Start Clearing !!!! \n ");
			for(int i=0;i<2;i++) {
				delete [] lclg_nonoverlap_transfrag[i];
				lclg_nonoverlap_transfrag[i] = new GPVec<CTransfrag>[20000];
			};
		}

		void Clear_no2gnode_unispg() {
		fprintf(stderr, "**** Start Clear_no2gnode_unispg !!!! \n ");
			for(int i=0;i<2;i++) {
				delete [] no2gnode_unispg[i];
				no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
		}

		void Clear_new_no2gnode_unispg() {
		// fprintf(stderr, "**** Start Clearing !!!! \n ");
			for(int i=0;i<2;i++) {
				delete [] new_no2gnode_unispg[i];
				new_no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
		}

		void Copy_new_no2gnode_unispg_2_no2gnode_unispg() {
		// fprintf(stderr, "**** Start Clearing !!!! \n ");

			fprintf(stderr, "**** Pre 'Copy_new_no2gnode_unispg_2_no2gnode_unispg' check !!!!\n");

			for(int s=0;s<2;s++) {
				
				// fprintf(stderr, ">> new_gidx[%d]: %d\n", s, new_gidx[s]);
				// for (int j=0; j<new_gidx[s]; j++) {
				fprintf(stderr, ">> 'Copy_new_no2gnode_unispg_2_no2gnode_unispg: 'new_no2gnode_unispg[%d][%d]: %d\n", s, 0, new_no2gnode_unispg[s]->Count());
				// GPVec<CGraphnodeUnispg> tmp = new GPVec(new_no2gnode_unispg[s][j]);
				no2gnode_unispg[s][current_gidx[s]] = new GPVec<CGraphnodeUnispg>(new_no2gnode_unispg[s][0]);
				// }
				// no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
			
			
			fprintf(stderr, "**** Post 'Copy_new_no2gnode_unispg_2_no2gnode_unispg' check !!!!\n");

			for(int s=0;s<2;s++) {
				
				//fprintf(stderr, ">> new_gidx[%d]: %d\n", s, new_gidx[s]);

				// for (int j=0; j<new_gidx[s]; j++) {
				fprintf(stderr, ">> 'Copy_new_no2gnode_unispg_2_no2gnode_unispg: 'no2gnode_unispg[%d][%d]: %d\n", s, current_gidx[s], no2gnode_unispg[s][current_gidx[s]].Count());
				for (int n=0; n<no2gnode_unispg[s][current_gidx[s]].Count(); n++) {
					fprintf(stderr, "\t>> no2gnode_unispg[%d][%d][%d]: %d\n", s, current_gidx[s], n, no2gnode_unispg[s][current_gidx[s]].Get(n)->nodeid);
				}
				// }
				// no2gnode_unispg[i] = new GPVec<CGraphnodeUnispg>[20000];
			};
		}

		void PrintGraphGp();
		GPVec<CGraphnodeUnispg>** get_no2gnodeGp ();
};

#endif