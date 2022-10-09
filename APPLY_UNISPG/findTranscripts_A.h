#ifndef __FINDTRANSCRIPTS_A_H__
#define __FINDTRANSCRIPTS_A_H__

#include "global_params_A.h"
#include "../unispg.h"
// #include "../helper.h"

// int find_transcripts(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand, BundleData* bdata);
int find_transcripts_APPLY_UNISPG(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand, GVec<CGuide>& guidetrf,GPVec<GffObj>& guides,GVec<int>& guidepred,BundleData* bdata,GVec<int>& trflong);

int guides_pushmaxflow_APPLY_UNISPG(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata);

void parse_trf_APPLY_UNISPG(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& usednode,float maxcov,GBitVec& prevpath);

float push_max_flow_APPLY_UNISPG(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,
		GVec<float>& nodeflux,GBitVec& pathpat, GIntHash<int> &gpos, bool &full);

float store_transcript_APPLY_UNISPG(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnodeUnispg>& no2gnode,int& geneno,bool& first,int strand,int gno,GIntHash<int>& gpos, bool& included,
		GBitVec& prevpath, bool full=false,BundleData *bdata=NULL, //float fragno, char* id=NULL) {
		   GffObj* t=NULL);

void update_guide_pred_APPLY_UNISPG(GList<CPrediction>& pred,int np, GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnodeUnispg>& no2gnode,int gno,bool update);

float push_guide_maxflow_APPLY_UNISPG(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,GBitVec& pathpat);

float guidepushflow_APPLY_UNISPG(int g,GVec<CGuide>& guidetrf,int gno,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnodeUnispg>& no2gnode,GVec<float>& nodeflux);

CTransfrag *find_guide_partial_pat_APPLY_UNISPG(GffObj *guide,GPVec<CGraphnodeUnispg>& no2gnode,int gno,int edgeno,GIntHash<int> &gpos,GVec<int>& olen,int &olensum);

int store_guide_transcript_APPLY_UNISPG(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnodeUnispg>& no2gnode,int& geneno,bool& first,int gno, GffObj* t,bool update);

bool back_to_source_fast_APPLY_UNISPG(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos);

bool fwd_to_sink_fast_APPLY_UNISPG(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos);

bool onpath_APPLY_UNISPG(GBitVec& trpattern,GVec<int>& trnode,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnodeUnispg>& no2gnode,int gno,
		GIntHash<int>& gpos);
#endif