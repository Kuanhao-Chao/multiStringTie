#ifndef __RLINK_C_H__
#define __RLINK_C_H__

#pragma once
#include "global_params_C.h"
#include "rlink_C_help.h"

CGraphnode *create_graphnode_unispg(int s, int g, uint start,uint end,int nodeno,CBundlenode *bundlenode,
		GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode);

GBitVec traverse_dfs_unispg(int s,int g,CGraphnode *node,CGraphnode *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos,int &lastgpos);

int create_graph_unispg(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode,
		GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GIntHash<int> **gpos,BundleData* bdata,
		int &edgeno,int &lastgpos,GArray<GEdge>& guideedge, int refend=0);

void process_transfrags_unispg(int s, int gno,int edgeno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no,
		GIntHash<int> &gpos,GVec<CGuide>& guidetrf,GList<CPrediction>& pred,GVec<int>& trflong);

int find_transcripts_unispg(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand,
		GVec<CGuide>& guidetrf,GPVec<GffObj>& guides,GVec<int>& guidepred,BundleData* bdata,GVec<int>& trflong);

CTreePat *construct_treepat_unispg(int gno, GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag);
#endif