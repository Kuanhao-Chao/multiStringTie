#include "rlink.h"
#include "rlink_multi.h"
#include "GBitVec.h"
#include <float.h>


#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif


extern FILE *c_out;         // file handle for the input transcripts that are fully covered by reads

extern GStr out_dir;
extern GStr outfname;
extern GStr plot_dir;
extern bool trim;
extern bool eonly;
extern bool nomulti;
extern bool viral;
extern bool mixedMode;
extern bool guided;

extern int allowed_nodes;
extern float isofrac;
extern bool isunitig;
extern bool longreads;
extern bool rawreads;
extern float mcov;
extern int mintranscriptlen; // minimum number for a transcript to be printed
extern uint junctionsupport; // anchor length for junction to be considered well supported <- consider shorter??
extern uint sserror;
extern int junctionthr; // number of reads needed to support a particular junction
extern float readthr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern float singlethr;     // read coverage per bundle bp to accept it; otherwise considered noise
extern uint bundledist;  // reads at what distance should be considered part of separate bundles
                        // <- this is not addressed everywhere, e.g. in infer_transcripts -> look into this

extern bool includesource;
extern bool geneabundance; // need to compute the gene abundance

extern float fpkm_thr;
extern float tpm_thr;
extern bool enableNames;
extern bool includecov;
extern bool retained_intron;

extern FILE* f_out;

/****************
 **  KH Adding 
****************/
extern FILE* uinigraph_out;
extern bool universal_splice_graph;
/****************
 **  END KH Adding 
****************/



extern GStr label;

inline int edge(int min, int max, int gno) {
	//return((gno-1)*min-min*(min-1)/2+max-min); // this should be changed if source to node edges are also stored
	return((gno-1)*(min+1)-min*(min-1)/2+max-min); // this includes source to node edges
}



CGraphnode *create_graphnode_multi(int s, int g, uint start,uint end,int nodeno,CBundlenode *bundlenode,
		GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"create_graphnode_multi[%d][%d]:%d-%d nodeno=%d\n",s,g,start,end,nodeno);
	}
	*/
/****************
 **  Original
 ****************/
	// CGraphnode* gnode=new CGraphnode(start,end,nodeno);
	// CGraphinfo ginfo(g,nodeno);
	// bundle2graph[s][bundlenode->bid].Add(ginfo);
	// no2gnode[s][g].Add(gnode);
/****************
 **  Original
 ****************/

/****************
 **  New
 ****************/
	CGraphnode* gnode=new CGraphnode(start,end,nodeno);
	CGraphinfo ginfo(g,nodeno);
	bundle2graph[s][bundlenode->bid].Add(ginfo);
	no2gnode[s][g].Add(gnode);
/****************
 **  New
 ****************/

	return(gnode);
}

CGraphnode *create_graphnode_multi_cov(int s, int g, uint start,uint end,int nodeno, float nodecov, CBundlenode *bundlenode, GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"create_graphnode_multi[%d][%d]:%d-%d nodeno=%d\n",s,g,start,end,nodeno);
	}
	*/
/****************
 **  Original
 ****************/
	// CGraphnode* gnode=new CGraphnode(start,end,nodeno);
	// CGraphinfo ginfo(g,nodeno);
	// bundle2graph[s][bundlenode->bid].Add(ginfo);
	// no2gnode[s][g].Add(gnode);
/****************
 **  Original
 ****************/

/****************
 **  New
 ****************/
	CGraphnode* gnode=new CGraphnode(start,end,nodeno, nodecov);
	CGraphinfo ginfo(g,nodeno);
	bundle2graph[s][bundlenode->bid].Add(ginfo);
	no2gnode[s][g].Add(gnode);
/****************
 **  New
 ****************/

	return(gnode);
}



CGraphnode *longtrim_multi(int s, int g, int refstart,int nodeend, int &nls, int &nle, bool &startcov, bool endcov, GVec<CPred> &lstart, GVec<CPred> &lend,
		CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>& futuretr, int& graphno, GVec<float>* bpcov,
		CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	while(nls<lstart.Count() && lstart[nls].predno<(int)graphnode->start) nls++;
	while(nle<lend.Count() && lend[nle].predno<(int)graphnode->start) nle++;
	while((nls<lstart.Count() && lstart[nls].predno<nodeend) || (nle<lend.Count() && lend[nle].predno<nodeend)){
		if(nle>=lend.Count() || (nls<lstart.Count() && lstart[nls].predno<=lend[nle].predno)) { // start comes first
			float tmpcov=0;
			if((startcov || lstart[nls].predno>(int)(graphnode->start+longintronanchor)) &&(endcov || lstart[nls].predno<nodeend+(int)longintronanchor)) { // start and ends can not be too close to a junction
				int startpos=lstart[nls].predno-refstart;
				int winstart=startpos-CHI_THR;
				if(winstart<0) winstart=0;
				int winend=startpos+CHI_THR-1;
				if(winend>=bpcov->Count()) winend=bpcov->Count()-1;
				/*tmpcov=(get_cov(1,startpos,winend,bpcov)-get_cov(2-2*s,startpos,winend,bpcov)-
						get_cov(1,winstart,startpos-1,bpcov)+get_cov(2-2*s,winstart,startpos-1,bpcov))/(DROP*CHI_THR);*/
				tmpcov=(get_cov_sign(2*s,startpos,winend,bpcov)-get_cov_sign(2*s,winstart,startpos-1,bpcov))/(DROP*CHI_THR);
			}
			if(tmpcov<=0 && lstart[nls].cov<0) tmpcov=ERROR_PERC; // to re-estimate later in process_transfrags
			if(tmpcov>0) {
				tmpcov+=trthr;
				uint prevend=graphnode->end;
				graphnode->end=lstart[nls].predno-1;
				CGraphnode *prevnode=graphnode;
				graphnode=create_graphnode_multi(s,g,lstart[nls].predno,prevend,graphno,bundlenode,bundle2graph,no2gnode);
				graphnode->hardstart=true;
/****************
 **  KH comment out
 ****************/
				graphno++;
				source->child.Add(graphnode->nodeid);  // this node is the child of source
				graphnode->parent.Add(source->nodeid); // this node has source as parent
				prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
				graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
/****************
 ** End of KH comment out
 ****************/
				float tmp=graphno-1;
				futuretr.cAdd(0.0);
				futuretr.Add(tmp);
				//tmp=lstart[nls].cov+trthr;
				futuretr.Add(tmpcov);
				//		graphnode->start,lstart[nls].cov,get_cov(2*s,startpos,startpos+CHI_THR-1,bpcov),get_cov(1,startpos,startpos+CHI_THR-1,bpcov),get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov),
				//		get_cov(2*s,startpos-CHI_THR,startpos-1,bpcov),get_cov(1,startpos-CHI_THR,startpos-1,bpcov),get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov));

				tmp=prevnode->nodeid;futuretr.Add(tmp);
				tmp=graphnode->nodeid;futuretr.Add(tmp);
				tmp=trthr;futuretr.Add(tmp);
				// COUNT 2 EDGES HERE
				edgeno+=2;
				startcov=false;
			}
			nls++;
		}
		else if(nls>=lstart.Count() || (nle<lend.Count() && lend[nle].predno<lstart[nls].predno)) { // end comes first
			float tmpcov=0;
			if((!startcov || lend[nle].predno>(int)(graphnode->start+longintronanchor)) &&(!endcov || lend[nle].predno<nodeend+(int)longintronanchor)) {
				int endpos=lend[nle].predno-refstart;
				int winstart=endpos-CHI_THR+1;
				if(winstart<0) winstart=0;
				int winend=endpos+CHI_THR;
				if(winend>=bpcov->Count()) winend=bpcov->Count()-1;
				/*tmpcov=(get_cov(1,winstart,endpos,bpcov)-get_cov(2-2*s,winstart,endpos,bpcov)-
						get_cov(1,endpos+1,winend,bpcov)+get_cov(2-2*s,endpos+1,winend,bpcov))/(DROP*CHI_THR);*/
				tmpcov=(get_cov_sign(2*s,winstart,endpos,bpcov)-get_cov_sign(2*s,endpos+1,winend,bpcov))/(DROP*CHI_THR);
			}
			if(tmpcov<=0 && lend[nle].cov<0) tmpcov=ERROR_PERC; // to re-estimate later in process_transfrags
			if(tmpcov>0) {
				tmpcov+=trthr;
				float tmp=graphno-1;
				uint prevend=graphnode->end;
				graphnode->end=lend[nle].predno;
				CGraphnode *prevnode=graphnode;
				graphnode->hardend=true;
				graphnode=create_graphnode_multi(s,g,lend[nle].predno+1,prevend,graphno,bundlenode,bundle2graph,no2gnode);
				graphno++;
/****************
 **  KH comment out
 ****************/
				prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
				graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
				sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
				// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
/****************
 ** End of KH comment out
 ****************/
				futuretr.Add(tmp);
				futuretr.cAdd(-1.0);
				//tmp=lend[nle].cov+trthr;
				futuretr.Add(tmpcov);
				tmp=prevnode->nodeid;futuretr.Add(tmp);
				tmp=graphnode->nodeid;futuretr.Add(tmp);
				tmp=trthr;futuretr.Add(tmp);
				// COUNT 2 EDGES HERE
				edgeno+=2;
				startcov=true;
			}
			nle++;
		}
	}
	return(graphnode);
}


CGraphnode *source2guide_multi(int s, int g, int refstart,uint newstart,uint newend, CGraphnode *graphnode,CGraphnode *source,
		GVec<float>* bpcov,GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode, int &edgeno) {

	if(graphnode->start+longintronanchor>newstart) { // newstart is very close to graphnode start
		for(int p=0;p<graphnode->parent.Count();p++) if(!graphnode->parent[p]) return(graphnode);
	}

	// compute maxabund
	float leftcov=0;
	float rightcov=0;

	//int os=2-2*s; // other strand

	if(!mergeMode) {
		if(newstart>graphnode->start) {
			uint gstart=graphnode->start;
			if(newstart-gstart > CHI_WIN) {
				gstart=newstart-CHI_WIN;
			}
			//leftcov=get_cov(1,gstart-refstart,newstart-1-refstart,bpcov)- get_cov(os,gstart-refstart,newstart-1-refstart,bpcov);
			leftcov=get_cov_sign(2*s,gstart-refstart,newstart-1-refstart,bpcov);
			leftcov/=newstart-gstart;
		}
		if(newstart<newend) {
			uint gend=newend;
			if(newend-newstart>=CHI_WIN) {
				gend=newstart+CHI_WIN-1;
			}
			//rightcov=get_cov(1,newstart-refstart,gend-refstart,bpcov)-get_cov(os,newstart-refstart,gend-refstart,bpcov);
			rightcov=get_cov_sign(2*s,newstart-refstart,gend-refstart,bpcov);
			rightcov/=gend-newstart+1;
		}
	}

	float maxabund=rightcov-leftcov;
	if(maxabund<trthr) maxabund=trthr;

	if(graphnode->start<=newstart-1) {
		uint prevend=graphnode->end;
		graphnode->end=newstart-1;
		CGraphnode *prevnode=graphnode;
		graphnode=create_graphnode_multi(s,g,newstart,prevend,graphno,bundlenode,bundle2graph,no2gnode);
		graphno++;
		float tmp=prevnode->nodeid;futuretr.Add(tmp);
		tmp=graphnode->nodeid;futuretr.Add(tmp);
		tmp=trthr;futuretr.Add(tmp);
/****************
 **  KH comment out
 ****************/
		prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
		graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
/****************
 ** End of KH comment out
 ****************/
	}
/****************
 **  KH comment out
 ****************/
	source->child.Add(graphnode->nodeid);  // this node is the child of source
	graphnode->parent.Add(source->nodeid); // this node has source as parent
/****************
 ** End of KH comment out
 ****************/
	float tmp=graphno-1;
	futuretr.cAdd(0.0);
	futuretr.Add(tmp);
	futuretr.Add(maxabund);
	// COUNT 1 EDGE HERE because the source to guide edge was already included in our count
	edgeno++;

	return(graphnode);

}


// cummulative version
CGraphnode *guide2sink_multi(int s, int g, int refstart,uint newstart,uint newend, CGraphnode *graphnode,CGraphnode *sink,
		GVec<float>* bpcov,GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode, int &edgeno) {

	// compute maxabund
	float leftcov=0;
	float rightcov=0;

	//int os=2-2*s;

	if(!mergeMode) {
		if(newstart>=graphnode->start) {
			uint gstart=graphnode->start;
			if(newstart-gstart >= CHI_WIN) {
				gstart=newstart-CHI_WIN+1;
			}
			//leftcov=get_cov(1,gstart-refstart,newstart-refstart,bpcov)-get_cov(os,gstart-refstart,newstart-refstart,bpcov);
			leftcov=get_cov_sign(2*s,gstart-refstart,newstart-refstart,bpcov);
			leftcov/=newstart-gstart+1;
		}
		if(newstart+1<newend) {
			uint gend=newend;
			if(newend-newstart>CHI_WIN) {
				gend=newstart+CHI_WIN;
			}
			//rightcov=get_cov(1,newstart+1-refstart,gend-refstart,bpcov)-get_cov(os,newstart+1-refstart,gend-refstart,bpcov);
			rightcov=get_cov_sign(2*s,newstart+1-refstart,gend-refstart,bpcov);
			rightcov/=gend-newstart;
		}
	}

	float maxabund=leftcov-rightcov;
	if(maxabund<trthr) maxabund=trthr;

	float tmp=graphno-1;
	uint prevend=graphnode->end;
	graphnode->end=newstart;
	CGraphnode *prevnode=graphnode;
	graphnode=create_graphnode_multi(s,g,newstart+1,prevend,graphno,bundlenode,bundle2graph,no2gnode);
	graphno++;
/****************
 **  KH comment out
 ****************/
	prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
	graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
	sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
/****************
 ** End of KH comment out
 ****************/
	futuretr.Add(tmp);
	futuretr.cAdd(-1.0);
	futuretr.Add(maxabund);
	tmp=prevnode->nodeid;futuretr.Add(tmp);
	tmp=graphnode->nodeid;futuretr.Add(tmp);
	tmp=trthr;futuretr.Add(tmp);
	// COUNT 1 EDGE HERE because the source to guide edge was already included in our count
	edgeno++;

	return(graphnode);

}


CGraphnode *trimnode_all_multi(int s, int g, int refstart,uint newend, CGraphnode *graphnode,CGraphnode *source, CGraphnode *sink, GVec<float>* bpcov,
		GVec<float>& futuretr, int& graphno,CBundlenode *bundlenode,GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode, int &edgeno) {

	GVec<CTrimPoint> trimpoint;
	find_all_trims(refstart,2*s,graphnode->start,newend,bpcov,trimpoint);
	for(int i=0;i<trimpoint.Count();i++) if(trimpoint[i].pos){
		if(trimpoint[i].start) { // source trim
			graphnode->end=trimpoint[i].pos-1;
			//fprintf(stderr,"Create source trim:%d-%d and %d-%d\n",graphnode->start,graphnode->end,trimpoint[i].pos,newend);
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode_multi(s,g,trimpoint[i].pos,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphnode->hardstart=true;
			graphno++;
/****************
 **  KH comment out
 ****************/
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
/****************
 ** End of KH comment out
 ****************/
			float tmp=graphno-1;
			futuretr.cAdd(0.0);
			futuretr.Add(tmp);
			float sourceabundance=trimpoint[i].abundance+trthr;futuretr.Add(sourceabundance);
			tmp=prevnode->nodeid;futuretr.Add(tmp);
			tmp=graphnode->nodeid;futuretr.Add(tmp);
			tmp=trthr;futuretr.Add(tmp);
			// COUNT 2 EDGES HERE
			edgeno+=2;
		}
		else { // this is sink
			graphnode->end=trimpoint[i].pos;
			CGraphnode *prevnode=graphnode;
			graphnode=create_graphnode_multi(s,g,trimpoint[i].pos+1,newend,graphno,bundlenode,bundle2graph,no2gnode);
			graphno++;
/****************
 **  KH comment out
 ****************/
			prevnode->child.Add(graphnode->nodeid); // this node is the child of previous node
			graphnode->parent.Add(prevnode->nodeid); // this node has as parent the previous node
			sink->parent.Add(prevnode->nodeid); // prevnode is the parent of sink
/****************
 ** End of KH comment out
 ****************/
			prevnode->hardend=true;
			// remember to create transfrag as well -> I don't know the gno yet, so I can not create it here
			float tmp=graphno-2;
			futuretr.Add(tmp);
			futuretr.cAdd(-1.0);
			float sinkabundance=trimpoint[i].abundance+trthr;futuretr.Add(sinkabundance);
			tmp=prevnode->nodeid;futuretr.Add(tmp);
			tmp=graphnode->nodeid;futuretr.Add(tmp);
			tmp=trthr;futuretr.Add(tmp);
			// COUNT 2 EDGES HERE
			edgeno+=2;
		}
	}

	return(graphnode);
}


GBitVec traverse_dfs_multi(int s,int g,CGraphnode *node,CGraphnode *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos,int &lastgpos){

	//fprintf(stderr,"Traverse node %d\n",node->nodeid);


	if(visit[node->nodeid]) {
		/*****************************
		 ** Step 1: Nodes have been visited
		 *****************************/
		node->parentpat = node->parentpat | parents;
		for(int n=0;n<gno;n++) {
			if(parents[n]) // add node's children to all parents of node
				no2gnode[s][g][n]->childpat = no2gnode[s][g][n]->childpat | node->childpat;
			else if(node->childpat[n])
				no2gnode[s][g][n]->parentpat = no2gnode[s][g][n]->parentpat | node->parentpat;
		}
	} else {
		/*****************************
		 ** Step 1: Nodes have not been visited
		 *****************************/
		node->childpat.resize(gno+edgeno);
		node->parentpat.resize(gno+edgeno);
		node->parentpat = node->parentpat | parents;
		visit[node->nodeid]=true;
		parents[node->nodeid]=1; // add the node to the parents

		/*****************************
		 ** node has source only as parent -> add transfrag from source to node.
		 *****************************/
		if(node->parent.Count()==1 && !node->parent[0]) {
			GBitVec trpat(gno+edgeno);
			trpat[0]=1;
			trpat[node->nodeid]=1;

			int key=edge(0,node->nodeid,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			else {
				gpos[s][g].Add(key,lastgpos);
				trpat[lastgpos]=1;
				lastgpos++;
			}

			GVec<int> nodes;
			nodes.cAdd(0);
			nodes.Add(node->nodeid);
			CTransfrag *tr=new CTransfrag(nodes,trpat,trthr);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add source transfrag[%d][%d]= %d and pattern",s,g,transfrag[s][g].Count());
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			*/
			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
			if(mixedMode) { // I need to add a long read as well
				CTransfrag *longtr=new CTransfrag(nodes,trpat,trthr);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			} else
			if(longreads) // || mixedMode)
				transfrag[s][g].Last()->longread=true;

		}

		/*****************************
		 ** node is not sink.
		 *****************************/
		int n=node->child.Count();
		if(node != sink && !n) {
			node->child.Add(sink->nodeid);  // add sink to the node's children
			sink->parent.Add(node->nodeid); // add node to sink's parents
			// create the transfrag that ends the node
			GBitVec trpat(gno+edgeno);
			trpat[node->nodeid]=1;
			trpat[gno-1]=1;

			int key=edge(node->nodeid,gno-1,gno);
			int *pos=gpos[s][g][key];
			if(pos!=NULL) trpat[*pos]=1;
			else {
				gpos[s][g].Add(key,lastgpos);
				trpat[lastgpos]=1;
				lastgpos++;
			}

			GVec<int> nodes;
			nodes.Add(node->nodeid);
			nodes.Add(sink->nodeid);
			CTransfrag *tr=new CTransfrag(nodes,trpat,trthr);

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add sink transfrag[%d][%d]= %d for nodeid=%d and pattern:",s,g,transfrag[s][g].Count(),node->nodeid);
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			*/

			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
			if(mixedMode) { // I need to add a long read as well
				CTransfrag *longtr=new CTransfrag(nodes,trpat,trthr);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			}
			else
			if(longreads)// || mixedMode)
				transfrag[s][g].Last()->longread=true;
			n++;
	    }
		/*
		fprintf(stderr,"Add %d children of node %d (%d-%d): ",n,node->nodeid,node->start,node->end);
		for(int i=0;i<n;i++) fprintf(stderr," %d",node->child[i]);
		fprintf(stderr,"\n");
		*/

		//edgeno+=n; // this will have to be deleted in the end; now I put it so that I can check equivalence with the one computed when creating the graph

		/*****************************
		 ** iterate through all children
		 *****************************/
	    for(int i=0; i< n; i++) {
	    	GBitVec childparents=parents;
	    	int min=node->nodeid; // nodeid is always smaller than child node ?
	    	int max=node->child[i];
	    	if(min>max) {
	    		max=node->nodeid; // nodeid is always smaller than child node ?
	    		min=node->child[i];
	    	}

			int key=edge(min,max,gno);
			int *pos=gpos[s][g][key];

			if(pos!=NULL) {
				childparents[*pos]=1; // add edge from node to child to the set of parents from child
				node->childpat[*pos]=1; // add edge from node to child to the set of node children
			}
			else {
				gpos[s][g].Add(key,lastgpos);
				childparents[lastgpos]=1;
				node->childpat[lastgpos]=1;
				lastgpos++;
			}

			//fprintf(stderr,"Call for child %d with id=%d\n",node->child[i],no2gnode[s][g][node->child[i]]->nodeid);
	    	node->childpat = node->childpat | traverse_dfs_multi(s,g,no2gnode[s][g][node->child[i]],sink,childparents,gno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	    }
	} // end else from if(visit[node->nodeid])

	GBitVec children = node->childpat;
	children[node->nodeid]=1;

	return(children);
}


int create_graph_multi(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode,
		GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GIntHash<int> **gpos,BundleData* bdata,
		int &edgeno,int &lastgpos,GArray<GEdge>& guideedge, UniSpliceGraphGp* uni_splice_graphGp, int refend=0){

	int uni_refstart = uni_splice_graphGp -> get_refstart();
	int uni_refend = uni_splice_graphGp -> get_refend();
	int* uni_gpSize = uni_splice_graphGp -> get_gpSize();
	GVec<int>* uni_graphnoGp = uni_splice_graphGp -> get_graphnoGp();
	GVec<int>* uni_edgenoGp = uni_splice_graphGp -> get_edgenoGp();  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
	GPVec<CGraphnode>** uni_no2gnodeGp = uni_splice_graphGp -> get_no2gnodeGp(); // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i

	// fprintf(stderr, "&&& g %d\n: ", g);
	// fprintf(stderr, "uni_refstart %d\n: ", uni_refstart);
	// fprintf(stderr, "uni_refend  %d\n: ", uni_refend);
	// fprintf(stderr, "uni_gpSize  %d\n: ", uni_gpSize[s]);
	// fprintf(stderr, "uni_graphnoGp  %d\n: ", uni_graphnoGp[s][g]);
	// fprintf(stderr, "uni_edgenoGp  %d\n: ", uni_edgenoGp[s][g]);

	/****************
	 **  KH Adding 
	****************/
	fprintf(stdout, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	fprintf(stdout, "&&&&&& Start 'create_graph_multi'\n");
	fprintf(stdout, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	/****************
	 **  END KH Adding 
	****************/

	GVec<float>* bpcov = bdata ? bdata->bpcov : NULL; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles

	// Traverse the graph
	// fprintf(stderr, ">> uni_gpSize[s]: %d\n", uni_gpSize[s]);
	// if(g <= uni_gpSize[s]) {
	// 	fprintf(stderr,"Digraph %d_%d_%d_%d {", uni_refstart, uni_refend, s, g);
	// 	// graphno[s][b]: number of nodes in graph.
	// 	if(uni_graphnoGp[s][g]) {
	// 		for(int nd=1;nd<uni_graphnoGp[s][g]-1;nd++)
	// 		fprintf(stderr,"%d[start=%d end=%d cov=%f];",nd,uni_no2gnodeGp[s][g][nd]->start,uni_no2gnodeGp[s][g][nd]->end,uni_no2gnodeGp[s][g][nd]->cov);

	// 	for(int nd=0;nd<uni_graphnoGp[s][g];nd++) {
	// 		// fprintf(stderr,"Node %d with parents:",i);
	// 		for(int c=0;c<uni_no2gnodeGp[s][g][nd]->child.Count();c++) {
	// 			fprintf(stderr,"%d->",nd);			
	// 			fprintf(stderr,"%d;",uni_no2gnodeGp[s][g][nd]->child[c]);
	// 		}
	// 	}
	// 	}
	// 	fprintf(stderr,"}\n");
	// }

/****************
 **  original
 ****************/
	// CGraphnode* source = uni_no2gnodeGp[g][0];
	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();
/****************
 **  original
 ****************/
	int njunctions=junction.Count();
	fprintf(stderr,"&&&&&&&&&&& Start graph[%d][%d] with %d edgeno and lastgpos=%d\n",s,g,edgeno,lastgpos);

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Junctions[%d][%d]: ",s,g);
		for(int i=0;i<njunctions;i++) fprintf(stderr," %d-%d:%d",junction[i]->start,junction[i]->end,junction[i]->strand);
		fprintf(stderr,"\n");
		fprintf(stderr,"eJunctions[%d][%d]: ",s,g);
		for(int i=0;i<njunctions;i++) fprintf(stderr," %d-%d:%d",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand);
		fprintf(stderr,"\n");

		// fprintf(stderr,"bundle2graph[%d][%d]: ",s,g);
		// for(int i=0;si<(*bundle2graph)->Count();i++) fprintf(stderr," %d-%d",(**bundle2graph)[i].ngraph,(**bundle2graph)[i].nodeno);
		// fprintf(stderr,"\n");
	}
	*/

	/*****************************
	 ** Step 1: Check start & end of bundle nodes (swap if necessary)
	 *****************************/
	int nge=0;
	bool processguide=false;
	CBundlenode *bundlenode=bnode[bundle->startnode];
	while(nge<guideedge.Count() && bundlenode!=NULL) {
		uint start=guideedge[nge].val;
		uint end=guideedge[nge].endval;
		if(start>end) Gswap(start,end);
		if(bundlenode->end<start) bundlenode=bundlenode->nextnode;
		else if(guideedge[nge].strand==s && bundlenode->start<=end) { // bundlenode->end>=start
			nge=0;
			processguide=true;
			break;
		}
		else nge++;
	}
	bundlenode=bnode[bundle->startnode];
	int njs=0; // index of sorted junction starts
	int nje=0; // index of sorted junction ends

	fprintf(stderr,"&&&&&&&&&&& process bundle %d-%d:%d bpcov_count=%d refstart=%d\n",bdata->start,bdata->end,s,bpcov->Count(),refstart);

/****************
 **  original
 ****************/
	int graphno=1; // number of nodes in graph
/****************
 **  original
 ****************/
	//GHash<GVec<int>* > ends; // keeps ids of all nodes ending at a certain position; OR ALL NODES THAT ARE LINKED BY JUNCTIONS TO A CERTAIN POSITION
    GIntHash< GVec<int>* > ends;
	GVec<float> futuretr; //future transfrags

	/*****************************
	 ** Step 2: I have a bunch of junctions at the start for which I need to create ends
	 *****************************/
	if(mergeMode) { // I have a bunch of junctions at the start for which I need to create ends

		while(njs<njunctions && !junction[njs]->start ) { // remember ends here for source node
			if((junction[njs]->strand+1) == 2*s) {
				//GStr je((int)junction[njs]->end);
				//GVec<int> *e=ends[je.chars()];
				GVec<int> *e=ends[junction[njs]->end];
				if(!e) {
					e = new GVec<int>();
					//ends.Add(je.chars(),e);
					ends.Add(junction[njs]->end, e);
				}
				e->cAdd(0);
			}
			njs++;
		}
	}

	/*****************************
	 ** Step 3: 'create_graphnode_multi' function
	 ** 	Process nodes in the bundle.
	 *****************************/
	int f=0; // feature index
	uint bundle_start=bundlenode->start;
	uint bundle_end=bnode[bundle->lastnodeid]->end;
	GHashMap<int, int> global2local_nodehash(false); //hash of pointers
	global2local_nodehash.Add(0, 0);


// I want to process bundles & compare the local to the global graph.
	int nd_global=1;


	vector<vector<int> > exonIntervals;
	vector<vector<int> > exonIntervals_unispg;
	vector<int> lcl_tmp;
	vector<int> unispg_tmp;
	// This is a local single bundle node. 

	int g_global = 0;
		// Get the first & last node of the universal splice graph
		// fprintf(stderr,"\n  &&& uni_no2gnodeGp[%d][%d][1] start: ", s, g_global);
		
		// The universal graph overlaps the local graph.
		// while (((unispg_start <= bundle_start) && (unispg_end >=bundle_start)) || ((unispg_start <= bundle_end) && (unispg_end >=bundle_end))) {
		// 	// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
		// 	// |(s)...............------(e)|-----    or   |(s).................(e)|----------   

		// }

	bool overlap = false;
	bool next_unispg = true;
	while (!overlap && next_unispg) {		
		// fprintf(stderr,"\n &&& Loop!! \n");
		if(g_global < uni_gpSize[s]) {
			// fprintf(stderr,"\n  &&& uni_graphnoGp[%d][%d]: %d \n", s, g_global, uni_graphnoGp[s][g_global]);
			if (uni_graphnoGp[s][g_global]>=3) {
				uint unispg_start = uni_no2gnodeGp[s][g_global][1]->start;
				uint unispg_end = uni_no2gnodeGp[s][g_global][uni_graphnoGp[s][g_global]-2]->end;
				// fprintf(stderr,"\n  &&& unispg_end: %d \n", unispg_end);
				// fprintf(stderr,"\n  &&& Inside while loop \n");
				fprintf(stderr,"\n  &&& unispg: %d - %d;  bundle: %d - %d \n", unispg_start, unispg_end, bundle_start, bundle_end);
				if (unispg_end < bundle_start) {
					// ----------   |(s).................(e)|
					fprintf(stderr,"\n  &&& Bundle: ----------   |(s).................(e)| \n");
					overlap = false;
					next_unispg = true;
				} else if (unispg_start < bundle_start && unispg_end >= bundle_start && unispg_end <= bundle_end) {
					// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
					fprintf(stderr,"\n  &&& Bundle: ----------|(s).................(e)|   or   -----|(s)-----............(e)| \n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start < bundle_start && unispg_end > bundle_end) {
					// -----|(s)------------(e)|--
					fprintf(stderr,"\n &&& Bundle: -----|(s)------------(e)|-- \n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start == bundle_start && unispg_end < bundle_end) {
					// |(s)----------.................(e)| 
					fprintf(stderr,"\n &&& Bundle: |(s)----------............(e)|\n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start > bundle_start && unispg_end < bundle_end) {
					// |(s)........----------........(e)|
					fprintf(stderr,"\n |(s)........----------........(e)| \n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start > bundle_start && unispg_end == bundle_end) {
					// |(s)............----------(e)|
					fprintf(stderr,"\n &&& Bundle: |(s)............----------(e)| \n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start == bundle_start && unispg_end == bundle_end) {
					// |(s)----------(e)|
					fprintf(stderr,"\n &&& Bundle: |(s)----------(e)| \n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start <= bundle_end && unispg_end > bundle_end) {
					// |(s)...............------(e)|-----    or   |(s).................(e)|----------   
					fprintf(stderr,"\n &&& Bundle: (s)...............------(e)|-----    or   |(s).................(e)|---------- \n");
					overlap = true;
					next_unispg = true;
				} else if (unispg_start > bundle_end) {
					// The node is outside the current bundle => This node belongs to the next bundlenode
					// |(s).................(e)|   ----------
					fprintf(stderr,"\n &&& Bundle: |(s).................(e)|   ---------- \n");
					overlap = false;
					next_unispg = false;
					break;
				} else {
					fprintf(stderr,"\n &&& Unknown area!!!! \n");
				}
				if (overlap) {
					fprintf(stderr,"Overlapping \n");
					for(int nd=1;nd<uni_graphnoGp[s][g_global]-1;nd++) {
						unispg_tmp.clear();
						unispg_tmp.push_back(uni_no2gnodeGp[s][g_global][nd]->start);
						unispg_tmp.push_back(uni_no2gnodeGp[s][g_global][nd]->end);
						exonIntervals_unispg.push_back(unispg_tmp);
					}
				}
				if (next_unispg) {
					fprintf(stderr,"Go to next universal splice graph \n");
					g_global += 1;
				}
			} else {
				g_global += 1;
			}
		} else {
			break;
		}
	}
// }



	while(bundlenode!=NULL) {
		lcl_tmp.clear();
		lcl_tmp.push_back(bundlenode->start);
		lcl_tmp.push_back(bundlenode->end);
		exonIntervals.push_back(lcl_tmp);
		// jhash.Add(ej, ejunction[j]);

		fprintf(stderr,"process bundlenode %d-%d:%d bpcov_count=%d refstart=%d\n",bundlenode->start,bundlenode->end,s,bpcov->Count(),refstart);




/****************
 **  This is the global graph
****************/
// Traverse the universal graph!!

// if(g < uni_gpSize[s]) {

// 	fprintf(stderr, "&&& s %d\n: ", s);
// 	fprintf(stderr, "&&& g %d\n: ", g);
// 	fprintf(stderr, "uni_refstart %d\n: ", uni_refstart);
// 	fprintf(stderr, "uni_refend  %d\n: ", uni_refend);
// 	fprintf(stderr, "uni_gpSize  %d\n: ", uni_gpSize[s]);
// 	fprintf(stderr, "uni_graphnoGp  %d\n: ", uni_graphnoGp[s][g]);
// 	fprintf(stderr, "uni_edgenoGp  %d\n: ", uni_edgenoGp[s][g]);

// 	fprintf(stderr,"nd_global: %d  Digraph %d_%d_%d_%d {", nd_global, uni_refstart, uni_refend, s, g);
// 	// graphno[s][b]: number of nodes in graph.

// // nd is the node index for the global graph.
// 	for(int nd=nd_global;nd<uni_graphnoGp[s][g]-1;nd++) {
// 		// fprintf(stderr,"nd: %d[start=%d end=%d cov=%f];",nd,uni_no2gnodeGp[s][g][nd]->start,uni_no2gnodeGp[s][g][nd]->end,uni_no2gnodeGp[s][g][nd]->cov);

// 		fprintf(stderr,"\n !!checker global splice graph node [start=%d end=%d cov=%f]; \n\tlocal bundle node [start=%d end=%d cov=%f];\n\n",uni_no2gnodeGp[s][g][nd]->start,uni_no2gnodeGp[s][g][nd]->end,uni_no2gnodeGp[s][g][nd]->cov, bundlenode->start, bundlenode->end);

// // The node belongs to the current bundlenode
// // bundlenode size >= global node size
// 		bool additional_node = false;
// 		bool global_local_overlap = false;
// 		if ((uni_no2gnodeGp[s][g][nd]->end < bundlenode->start)) {
// 			additional_node = true;
// 			// ----------   |(s).................(e)|
// 			fprintf(stderr,"\n  &&& Bundle: ----------   |(s).................(e)| \n");
// 			// continue;
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start < bundlenode->start && uni_no2gnodeGp[s][g][nd]->end >= bundlenode->start && uni_no2gnodeGp[s][g][nd]->end <= bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
// 			fprintf(stderr,"\n  &&& Bundle: ----------|(s).................(e)|   or   -----|(s)-----............(e)| \n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start < bundlenode->start && uni_no2gnodeGp[s][g][nd]->end > bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// -----|(s)------------(e)|--
// 			fprintf(stderr,"\n &&& Bundle: -----|(s)------------(e)|-- \n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start == bundlenode->start && uni_no2gnodeGp[s][g][nd]->end < bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// |(s)----------.................(e)| 
// 			fprintf(stderr,"\n &&& Bundle: |(s)----------............(e)|\n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start > bundlenode->start && uni_no2gnodeGp[s][g][nd]->end < bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// |(s)........----------........(e)|
// 			fprintf(stderr,"\n |(s)........----------........(e)| \n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start > bundlenode->start && uni_no2gnodeGp[s][g][nd]->end == bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// |(s)............----------(e)|
// 			fprintf(stderr,"\n &&& Bundle: |(s)............----------(e)| \n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start == bundlenode->start && uni_no2gnodeGp[s][g][nd]->end == bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// |(s)----------(e)|
// 			fprintf(stderr,"\n &&& Bundle: |(s)----------(e)| \n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start <= bundlenode->end && uni_no2gnodeGp[s][g][nd]->end > bundlenode->end) {
// 			additional_node = true;
// 			global_local_overlap = true;
// 			// |(s)...............------(e)|-----    or   |(s).................(e)|----------   
// 			fprintf(stderr,"\n &&& Bundle: (s)...............------(e)|-----    or   |(s).................(e)|---------- \n");
// 			nd_global = nd+1;
// 		} else if (uni_no2gnodeGp[s][g][nd]->start > bundlenode->end) {
// 			// The node is outside the current bundle => This node belongs to the next bundlenode
// 			// |(s).................(e)|   ----------
// 			fprintf(stderr,"\n &&& Bundle: |(s).................(e)|   ---------- \n");
// 			// I should create another bundle node!!!
// 			// bundlenode=bundlenode->nextnode; // advance to next bundlenode
// 			break;
// 		}
// 		if (additional_node) {
// 			unispg_tmp.clear();
// 			unispg_tmp.push_back(uni_no2gnodeGp[s][g][nd]->start);
// 			unispg_tmp.push_back(uni_no2gnodeGp[s][g][nd]->end);
// 			exonIntervals_unispg.push_back(unispg_tmp);
// 		// 	CGraphnode *graphnode=create_graphnode_multi_cov(s,g,uni_no2gnodeGp[s][g][nd]->start,uni_no2gnodeGp[s][g][nd]->end,graphno,uni_no2gnodeGp[s][g][nd]->cov+1,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
// 		// 	fprintf(stderr, "\nAdding !!! ^^^ global: %d, local: %d.  \n", nd, graphno);
// 		// 	global2local_nodehash.Add(nd, graphno);
// 		// 	// I need to create a mapping from old node id to new node id!!!
// 		// 	// Adding edges here
// 		// 	for(int p=0;p<uni_no2gnodeGp[s][g][nd]->parent.Count();p++) {
// 		// 		fprintf(stderr,"%d->%d;", uni_no2gnodeGp[s][g][nd]->parent[p], nd);
// 		// 		const int* local_node_idx=global2local_nodehash[uni_no2gnodeGp[s][g][nd]->parent[p]];
// 		// 		if (local_node_idx) {
// 		// 			fprintf(stderr, "\nRetrieving!! ^^^ global: %d, local: %d.  \n", uni_no2gnodeGp[s][g][nd]->parent[p], *local_node_idx);
// 		// 			// int local_node_idx = global2local_nodehash[uni_no2gnodeGp[s][g][nd]->parent[p]];
// 		// 			// In a matched graph, now I need to add edges.
// 		// 			// CGraphnode *node_p=no2gnode[s][g][uni_no2gnodeGp[s][g][nd]->parent[p]];



// 		// 			graphnode -> parent.Add(uni_no2gnodeGp[s][g][nd]->parent[p]);

// 		// 			CGraphnode *node_p=no2gnode[s][g][uni_no2gnodeGp[s][g][nd]->parent[p]];
// 		// 			node_p -> child.Add(graphnode->nodeid);

// 		// 			// COUNT EDGE HERE
// 		// 			edgeno++;
// 		// 			// fprintf(stderr,"1 Edge %d-%d, edgeno=%d\n",node_p->nodeid,graphnode->nodeid,edgeno);
// 		// 		}
// 		// 	}
// 		// 	graphno++;

// 		// 	// This is the sink of the universal splice graph.
// 		// 	// CGraphnode *unispg_sink = uni_no2gnodeGp[s][g][uni_graphnoGp[s][g]-1];
			
// 		// 	for(int c=0;c<uni_no2gnodeGp[s][g][nd]->child.Count();c++) {
// 		// 		if (uni_no2gnodeGp[s][g][nd]->child[c] == uni_graphnoGp[s][g]-1) {
// 		// 			fprintf(stderr, "This node connect to the sink!!!\n");
// 		// 			// graphnode -> child.Add(sink->nodeid);
// 		// 			sink->parent.Add(graphnode->nodeid);
// 		// 			edgeno++;
// 		// 		}

// 		// 	}
// 		}  // end of => if (additional_node) {	
// 	}




// 	// for(int nd=0;nd<uni_graphnoGp[s][g];nd++) {
// 	// 	// fprintf(stderr,"Node %d with parents:",i);
// 	// 	for(int c=0;c<uni_no2gnodeGp[s][g][nd]->child.Count();c++) {
// 	// 		fprintf(stderr,"%d->",nd);			
// 	// 		fprintf(stderr,"%d;",uni_no2gnodeGp[s][g][nd]->child[c]);
// 	// 	}
// 	// }

// 	// 	for(int nd=0;nd<uni_graphnoGp[s][g];nd++) {
// 	// 		// fprintf(stderr,"Node %d with parents:",i);
// 	// 		for(int c=0;c<uni_no2gnodeGp[s][g][nd]->child.Count();c++) {
// 	// 			fprintf(stderr,"%d->",nd);			
// 	// 			fprintf(stderr,"%d;",uni_no2gnodeGp[s][g][nd]->child[c]);
// 	// 		}
// 	// 	}
// 	fprintf(stderr,"}\n");
// }

// /****************
//  **  This is the global graph (End)
// ****************/


// 					// We need to consider the global graph and the local junctions distribution.
// // We need to massage the node!
// 					uint currentstart=bundlenode->start; // current start is bundlenode's start
// 					uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now

// 					int end=0;

// 					while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
// 						if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
// 							end=1;
// 						}
// 						nje++;
// 					}

// 					GVec<CPred> lstart; // CPred: prediction point class
// 					GVec<CPred> lend;
// 					int fs=-1; // first start feature index in lstart
// 					int fe=-1; // first end feature index in lend

// 					// see if I need to adjust the start to ignore little hanging pieces that make no sense
// 					if(!end) {
// 						while(nje<njunctions && ejunction[nje]->strand+1!=2*s) nje++; // skip all junctions that are not on the same strand
// 						if(!mergeMode && (nje<njunctions && ejunction[nje]->end - currentstart < junctionsupport) &&
// 								(fs<0 || (uint)lstart[fs].predno>=ejunction[nje]->end) &&  // I do not want to miss any hard starts/ends
// 								(fe<0 || (uint)lend[fe].predno>=ejunction[nje]->end)) { // there is a junction ending soon here
// 							float covleft=get_cov(1,currentstart-refstart,ejunction[nje]->end-1-refstart,bpcov);
// 							float covright=get_cov(1,ejunction[nje]->end-refstart,2*ejunction[nje]->end - currentstart-1-refstart,bpcov);
// 							if(covleft<covright*(1-ERROR_PERC)) { // adjust start here if needed
// 								currentstart=ejunction[nje]->end;
// 								// I have to check ending junctions here again
// 								while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
// 									if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
// 										end=1;
// 									}
// 									nje++;
// 								}
// 							}
// 						}
// 					}

					
// 			/****************
// 			 **  Creating node based on the universal splice graph
// 			****************/
// 					// CGraphnode *graphnode=create_graphnode_multi_cov(s,g,uni_no2gnodeGp[s][g][nd]->start,uni_no2gnodeGp[s][g][nd]->end,graphno,uni_no2gnodeGp[s][g][nd]->cov+1,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
// 					// fprintf(stderr, "\nAdding !!! ^^^ global: %d, local: %d.  \n", nd, graphno);
// 					// global2local_nodehash.Add(nd, graphno);
// 					// // I need to create a mapping from old node id to new node id!!!
// 					// // Adding edges here
// 					// for(int p=0;p<uni_no2gnodeGp[s][g][nd]->parent.Count();p++) {
// 					// 	fprintf(stderr,"%d->%d;", uni_no2gnodeGp[s][g][nd]->parent[p], nd);
// 					// 	const int* local_node_idx=global2local_nodehash[uni_no2gnodeGp[s][g][nd]->parent[p]];
// 					// 	if (local_node_idx) {
// 					// 		fprintf(stderr, "\nRetrieving!! ^^^ global: %d, local: %d.  \n", uni_no2gnodeGp[s][g][nd]->parent[p], *local_node_idx);
// 					// 		// int local_node_idx = global2local_nodehash[uni_no2gnodeGp[s][g][nd]->parent[p]];
// 					// 		// In a matched graph, now I need to add edges.
// 					// 		// CGraphnode *node_p=no2gnode[s][g][uni_no2gnodeGp[s][g][nd]->parent[p]];
// 					// 		graphnode -> parent.Add(uni_no2gnodeGp[s][g][nd]->parent[p]);
// 					// 		CGraphnode *node_p=no2gnode[s][g][uni_no2gnodeGp[s][g][nd]->parent[p]];
// 					// 		node_p -> child.Add(graphnode->nodeid);
// 					// 		// COUNT EDGE HERE
// 					// 		edgeno++;
// 					// 		// fprintf(stderr,"1 Edge %d-%d, edgeno=%d\n",node_p->nodeid,graphnode->nodeid,edgeno);
// 					// 	}
// 					// }
// 					// graphno++;
// 					// nd_global = nd+1;
// 					// fprintf(stderr,"create graph 1\n");
// 			/****************
// 			 **  Creating node based on the universal splice graph
// 			****************/

// 			/****************
// 			 **  original
// 			****************/
// 					CGraphnode *graphnode=create_graphnode_multi(s,g,currentstart,endbundle,graphno,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
// 					graphno++;
// 			/****************
// 			 **  original
// 			****************/

// 			// KH: For multi-samples, I should just skip this.
// 			/****************
// 			 **  original
// 			****************/
// 					if(end) { // I might have nodes finishing here; but I have a junction finishing here for sure
// 						//GStr cs((int)currentstart);
// 						//GVec<int> *e=ends[cs.chars()]; // HOW CAN I HAVE MORE THAN ONE NODE FINISHING HERE???; because this keeps all nodes that are linked by junctions here
// 						GVec<int> *e=ends[currentstart];
// 						if(e) {
// 							for(int i=0;i<e->Count();i++) {
// 								CGraphnode *node=no2gnode[s][g][e->Get(i)];
// 								node->child.Add(graphnode->nodeid);  // this node is the child of previous node
// 								graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
// 								// COUNT EDGE HERE
// 								edgeno++;
// 								fprintf(stderr,"1 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
// 							}
// 						}
// 						else { // I haven't seen nodes before that finish here (maybe due to error correction?) => link to source
// 					    	source->child.Add(graphnode->nodeid);  // this node is the child of source
// 					    	graphnode->parent.Add(source->nodeid); // this node has source as parent
// 					    	// COUNT EDGE HERE
// 					    	edgeno++;
// 					    	fprintf(stderr,"2 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
// 						}
// 					}
// 					else { // this node comes from source directly
// 						source->child.Add(graphnode->nodeid);  // this node is the child of source
// 						graphnode->parent.Add(source->nodeid); // this node has source as parent
// 						// COUNT EDGE HERE
// 						edgeno++;
// 						fprintf(stderr,"3 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
// 					}
// 			/****************
// 			 **  original
// 			****************/

// 					bool completed=false;

// 					bool dropcov=false; // false(0) means start of bundle or junction end (raise in coverage); true(1) means junction start (drop in coverage)
// 					int nls=0; // index in longstart
// 					int nle=0; // index in longend


// 					do {
// 						while(nje<njunctions && (((int)ejunction[nje]->strand+1) != 2*s)) nje++; // skip junctions that don't have the same strand
// 						while(njs<njunctions && ((((int)junction[njs]->strand+1) != 2*s) || (junction[njs]->start<currentstart))) njs++; // junctions that start before the current graphnode and I haven't seen them before are part of a different bundle


// 						int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome; if they have same position then process ending first
// 						if((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
// 							if(njs<njunctions && (junction[njs]->start<=endbundle) && junction[njs]->end>bundle_end) njs++;
// 							else {
// 								if(nje<njunctions) { // there are still junctions endings
// 									if(njs<njunctions) { // there are still junctions starting
// 										minjunction = junction[njs]->start >= ejunction[nje]->end ? 1 : 0; // one of them is clearly before the endbundle from the initial if
// 									}
// 									else minjunction = 1;
// 								}
// 								else minjunction = 0;
// 							}
// 						}

// 						// fprintf(stderr,"minjunction=%d\n",minjunction);
// 						if(nje<njunctions) fprintf(stderr,"Found junction:%d-%d(%d)\n",ejunction[nje]->start,ejunction[nje]->end,ejunction[nje]->strand);


// 			/****************
// 			**  original
// 			****************/
// 						if(minjunction == 0 ) { // found a start junction here

// 		// add guide starts/ends first
// 		if(processguide) {
// 			while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
// 			if(nge<guideedge.Count()) {

// 				while(true) {

// 					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

// 					if(nge>=guideedge.Count() || guideedge[nge].val>=junction[njs]->start) break;

// 					uint gstart=guideedge[nge].val;
// 					uint gend=junction[njs]->start;
// 					bool sourceguide=false;
// 					if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
// 					nge++;
// 					if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
// 					else if(guideedge[nge-1].endval<currentstart) continue;

// 					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
// 					if(nge<guideedge.Count() && guideedge[nge].val<junction[njs]->start) gend=guideedge[nge].val;

// 					// I need to check there is no other trimming needed due to drops from longreads
// 					if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,gstart,nls,nle,dropcov,!sourceguide,lstart,lend,
// 							graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

// 					if(sourceguide)	{
// 						graphnode=source2guide(s,g,refstart,gstart,gend,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
// 						dropcov=false;
// 					}
// 					else {
// 						graphnode=guide2sink(s,g,refstart,gstart,gend,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
// 						dropcov=true;
// 					}

// 				}
// 			}
// 		}
// 		// if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes
// 		else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,junction[njs]->start,nls,nle,dropcov,true,lstart,lend,
// 					graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
// 		if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


// 							dropcov=true;
// 							// if no trimming required just set the end of the node
// 							graphnode->end=junction[njs]->start; // set the end of current graphnode to here; introduce smaller nodes if trimming is activated
// 							uint pos=junction[njs]->start;
// 							while(njs<njunctions && junction[njs]->start==pos ) { // remember ends here
// 			/****************
// 			**  KH comment out
// 			****************/
// 								if((junction[njs]->strand+1) == 2*s) {
// 									//seenjunc++;
// 									if(mergeMode && (int)junction[njs]->end==refend) { // this node goes straight to sink
// 										sink->parent.Add(graphnode->nodeid); // graphnode is the parent of sink: check to see if I have a conflict with this
// 										edgeno++; // count edge here
// 									}
// 									else {
// 										//GStr je((int)junction[njs]->end);
// 										//GVec<int> *e=ends[je.chars()];
// 										GVec<int> *e=ends[junction[njs]->end];
// 										if(!e) {
// 											e = new GVec<int>();
// 											//ends.Add(je.chars(),e);
// 											ends.Add(junction[njs]->end, e);
// 										}
// 										e->Add(graphnode->nodeid);
// 									}
// 								}
// 			/****************
// 			** End of KH comment out
// 			****************/
// 								njs++;
// 							}

// 			/****************
// 			**  KH comment out
// 			****************/
// 				    		if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
// 				    			// see if I should just skip node
// 				    			if(endbundle-pos<junctionsupport) {
// 				    				while(njs<njunctions && junction[njs]->strand+1 != 2*s) njs++;
// 				    				if(!mergeMode && (njs>=njunctions || junction[njs]->start > endbundle) && (nje>=njunctions || ejunction[nje]->end > endbundle)) { // there are no more junctions starting within this bundle
// 				    					float covleft=get_cov(1,2*pos-endbundle+1-refstart,pos-refstart,bpcov);
// 				    					float covright=get_cov(1,pos+1-refstart,endbundle-refstart,bpcov);
// 				    					if(covright<covleft*(1-ERROR_PERC)) { // adjust start here if needed
// 				    						completed=true;
// 				    					}
// 				    				}
// 				    			}

// 				    			if(!completed) {
// 				    				//fprintf(stderr,"create graph 2\n");
// 				    				// CGraphnode *nextnode = create_graphnode_multi_cov(s,g,pos+1,endbundle,graphno,1,bundlenode,bundle2graph,no2gnode);
// 	    							CGraphnode *nextnode = create_graphnode_multi(s,g,pos+1,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
// 				    				graphno++;
// 				    				graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
// 				    				nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
// 				    				// COUNT EDGE HERE
// 				    				edgeno++;
// 				    				fprintf(stderr,"4 Edge %d-%d, edgeno=%d nextnode: %u-%u pos=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno,nextnode->start,nextnode->end,pos);
// 				    				graphnode=nextnode;
// 				    			}
// 				    		}
// 				    		else completed=true;
// 			/****************
// 			** End of KH comment out
// 			****************/
// 						}
// 						else if(minjunction == 1) { // found a junction end here

// 							uint pos=ejunction[nje]->end;
// 							while(nje<njunctions && ejunction[nje]->end==pos) { // read all junction ends at the current start
// 								nje++;
// 							}

// 			/****************
// 			**  KH comment out
// 			****************/
// 							if(graphnode->start<pos) { // last created node starts before the position of the new node I want to create




// 		// add guide starts/ends first
// 		if(processguide) {
// 			while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
// 			if(nge<guideedge.Count()) {

// 				while(true) {

// 					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

// 					if(nge>=guideedge.Count() || guideedge[nge].val>=pos-1) break;

// 					uint start=guideedge[nge].val;
// 					uint end=pos-1;
// 					bool sourceguide=false;
// 					if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
// 					nge++;
// 					if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
// 					else if(guideedge[nge-1].endval<currentstart) continue;

// 					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
// 					if(nge<guideedge.Count() && guideedge[nge].val<pos-1) end=guideedge[nge].val;

// 					// I need to check there is no other trimming needed due to drops from longreads
// 					if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,start,nls,nle,dropcov,!sourceguide,lstart,lend,
// 							graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

// 					if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
// 					else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

// 				}
// 			}
// 		}
// 		//if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes
// 		else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,pos-1,nls,nle,dropcov,false,lstart,lend,
// 					graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
// 		if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes





// 								graphnode->end=pos-1; // set end of current graphnode here
// 								dropcov=false;
// 								// fprintf(stderr,"create graph 3\n");
// 								// CGraphnode *nextnode = create_graphnode_multi_cov(s,g,pos,endbundle,graphno,1,bundlenode,bundle2graph,no2gnode);
// 	    						CGraphnode *nextnode = create_graphnode_multi(s,g,pos,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
// 								graphno++;
// 								graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
// 								nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode

// 								// COUNT EDGE HERE
// 								edgeno++;
// 								fprintf(stderr,"5 Edge %d-%d, edgeno=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno);

// 								graphnode=nextnode;
// 							}
// 			/****************
// 			** End of KH comment out
// 			****************/	

// 			/****************
// 			**  KH comment out
// 			****************/
// 							// GStr spos((int)pos);
// 							// GVec<int> *e=ends[spos.chars()]; // WHY DOESN'T THIS REPEAT THE SAME THING IN CASE THE START HASN'T BEEN ADJUSTED? because nje is bigger now than the ones that end at the currentstart

// 							GVec<int> *e=ends[pos];
// 							if(e) for(int i=0;i<e->Count();i++) {
// 								CGraphnode *node=no2gnode[s][g][e->Get(i)];
// 								node->child.Add(graphnode->nodeid);  // this node is the child of previous node
// 								graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
// 								// COUNT EDGE HERE
// 								edgeno++;
// 								fprintf(stderr,"6 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
// 							}
// 			/****************
// 			** End of KH comment out
// 			****************/	
// 						}
// 			/****************
// 			**  original
// 			****************/


// 					} while((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle)));


// 			/****************
// 			**  KH comment out
// 			****************/
// 					if(!completed) { // I did not finish node --> this will be an ending node





// 		// add guide starts/ends first
// 		if(processguide) {
// 			while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
// 			if(nge<guideedge.Count()) {

// 				while(true) {

// 					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

// 					if(nge>=guideedge.Count() || guideedge[nge].val>=endbundle) break;

// 					uint start=guideedge[nge].val;
// 					uint end=endbundle;
// 					bool sourceguide=false;
// 					if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
// 					nge++;
// 					if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
// 					else if(guideedge[nge-1].endval<currentstart) continue;

// 					while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
// 					if(nge<guideedge.Count() && guideedge[nge].val<endbundle) end=guideedge[nge].val;

// 					// I need to check there is no other trimming needed due to drops from longreads
// 					if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,start,nls,nle,dropcov,!sourceguide,lstart,lend,
// 							graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

// 					if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
// 					else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

// 				}
// 			}
// 		}
// 		// if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes
// 		else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,endbundle,nls,nle,dropcov,true,lstart,lend,
// 					graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
// 		if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes





// 						graphnode->end=endbundle;
// 						// COUNT EDGE HERE (this is an edge to sink)
// 						edgeno++;
// 						fprintf(stderr,"7 Edge to sink from %d, edgeno=%d\n",graphnode->nodeid,edgeno);
// 					}
// 			/****************
// 			** End of KH comment out
// 			****************/	

// /****************
//  **  This is the start of global graph (Lower part)
// ****************/
				
// 		// }
// /****************
//  **  This is the end of the global graph (Lower part)
// ****************/

	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)


	if (!exonIntervals.empty() && !exonIntervals_unispg.empty()) {
	// if (!exonIntervals.empty()) {
		for (int i = 0; i < exonIntervals.size(); i++) {
			fprintf(stderr, "exonIntervals[%d][%d]: %d\n", i, 0,exonIntervals[i][0]);
			fprintf(stderr, "exonIntervals[%d][%d]: %d\n", i, 1,exonIntervals[i][1]);
			float span = static_cast<float>(exonIntervals[i][1] - exonIntervals[i][0]);
			fprintf(stderr, "span: %f\n", span);
		}
		if (exonIntervals_unispg.empty()) {
			// draw(exonIntervals, exonIntervals, "./plot_"+to_string(uni_refstart)+"_"+to_string(uni_refend)+"_"+to_string(s)+"_"+to_string(g)+".png");
		} else {
			for (int i = 0; i < exonIntervals_unispg.size(); i++) {
				fprintf(stderr, "exonIntervals_unispg[%d][%d]: %d\n", i, 0,exonIntervals_unispg[i][0]);
				fprintf(stderr, "exonIntervals_unispg[%d][%d]: %d\n", i, 1,exonIntervals_unispg[i][1]);
				float span = static_cast<float>(exonIntervals_unispg[i][1] - exonIntervals_unispg[i][0]);
				fprintf(stderr, "span: %f\n", span);
			}

			if (exonIntervals.size() >= 4 && exonIntervals_unispg.size() >= 4) {

				fprintf(stderr, "plot_dir: %s\n", plot_dir.chars());
				

				draw(exonIntervals, exonIntervals_unispg, string(plot_dir.chars())+"/plot_"+to_string(uni_refstart)+"_"+to_string(uni_refend)+"_"+to_string(s)+"_"+to_string(g)+".png");
			}
		}	
	}
    


// Traverse the graph
	if(no2gnode[s][g].Count()) {

		fprintf(stderr,"Traversing the created graph!!!\n");
		fprintf(stderr,"Digraph %d_%d_%d_%d {", bdata->start,bdata->end, s, g);
		// graphno[s][b]: number of nodes in graph.
		if(no2gnode[s][g].Count()) {
			for(int nd=1;nd<no2gnode[s][g].Count()-1;nd++)
			fprintf(stderr,"%d[start=%d end=%d cov=%f];",nd,no2gnode[s][g][nd]->start,no2gnode[s][g][nd]->end,no2gnode[s][g][nd]->cov);

		for(int nd=0;nd<no2gnode[s][g].Count();nd++) {
			// fprintf(stderr,"Node %d with parents:",i);
			for(int c=0;c<no2gnode[s][g][nd]->child.Count();c++) {
				fprintf(stderr,"%d->",nd);			
				fprintf(stderr,"%d;",no2gnode[s][g][nd]->child[c]);
			}
		}
		}
		fprintf(stderr,"}\n");
	}

	// fprintf(stderr,"graphno=%d\n",graphno);

	/*****************************
	 ** Step 4: 'get_cov_sign' function
	 ** 	add source/sink links for very high coverage drops
	 *****************************/
	if(!mergeMode) for(int i=1;i<graphno;i++) {
		float icov=0;
		if(i>1 && no2gnode[s][g][i]->parent[0]) { // node i has parents, and does not come from source => might need to be linked to source
			/*icov=(get_cov(1,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)-
					get_cov(2-2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov))/no2gnode[s][g][i]->len();*/
			icov=get_cov_sign(2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)/no2gnode[s][g][i]->len();
			float parcov=0; // all parents' coverage
			bool addsource=true;
			for(int j=0;j<no2gnode[s][g][i]->parent.Count();j++) {
				if(!no2gnode[s][g][i]->parent[j]) { // parent is source
					addsource=false;break;
				}
				int p=no2gnode[s][g][i]->parent[j];
				/*parcov+=(get_cov(1,no2gnode[s][g][p]->start-refstart,no2gnode[s][g][p]->end-refstart,bpcov)-
						get_cov(2-2*s,no2gnode[s][g][p]->start-refstart,no2gnode[s][g][p]->end-refstart,bpcov))/no2gnode[s][g][p]->len();*/
				parcov+=get_cov_sign(2*s,no2gnode[s][g][p]->start-refstart,no2gnode[s][g][p]->end-refstart,bpcov)/no2gnode[s][g][p]->len();
			}
			if(addsource && parcov<icov*ERROR_PERC*DROP) {
				no2gnode[s][g][i]->parent.Insert(0,source->nodeid);
				source->child.Add(i);
				futuretr.cAdd(0.0);
				float tmp=i;
				futuretr.Add(tmp);
				tmp=(icov-parcov)/DROP;
				futuretr.Add(tmp);
				// fprintf(stderr,"Add link to source from node %d with abundance %f\n",i,tmp);
				edgeno++;
			}
		}

		if(i<graphno-1) {
			if(no2gnode[s][g][i]->child.Count()) { // might need to be linked to sink
				bool addsink=true;
				for(int j=0;j<sink->parent.Count();j++) {
					if(sink->parent[j]==i) {
						addsink=false;
						break;
					}
				}
				if(addsink) {
					/*if(!icov) icov=(get_cov(1,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)-
							get_cov(2-2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov))/no2gnode[s][g][i]->len();*/
					if(!icov) icov=get_cov_sign(2*s,no2gnode[s][g][i]->start-refstart,no2gnode[s][g][i]->end-refstart,bpcov)/no2gnode[s][g][i]->len();
					float chcov=0; // all children coverages
					for(int j=0;j<no2gnode[s][g][i]->child.Count();j++) {
						int c=no2gnode[s][g][i]->child[j];
						/*chcov+=(get_cov(1,no2gnode[s][g][c]->start-refstart,no2gnode[s][g][c]->end-refstart,bpcov)-
								get_cov(2-2*s,no2gnode[s][g][c]->start-refstart,no2gnode[s][g][c]->end-refstart,bpcov))/no2gnode[s][g][c]->len();*/
						chcov+=get_cov_sign(2*s,no2gnode[s][g][c]->start-refstart,no2gnode[s][g][c]->end-refstart,bpcov)/no2gnode[s][g][c]->len();
					}
					if(chcov<icov*ERROR_PERC*DROP) {
						sink->parent.Add(i); // prevnode is the parent of sink
						float tmp=i;
						futuretr.Add(tmp);
						futuretr.cAdd(-1.0);
						tmp=(icov-chcov)/DROP;
						futuretr.Add(tmp);
						edgeno++;
					}
				}
			}
			else edgeno++; // node has no children so it might get linked to sink in traverse_dfs_multi function
		}
	}

	/*****************************
	 ** Step 5: 'prune_graph_nodes' function
	 ** 	here I know graphno => I can see if it's too big
	 *****************************/
	if(!mergeMode && graphno>allowed_nodes) { // TODO: define allowed_nodes as a default in stringtie.cpp that varies with the memory
		graphno=prune_graph_nodes(graphno,s,g,bundle2graph,bnode.Count(),no2gnode,junction,edgeno,futuretr,sink);
	}

/****************
 **  `no2gnode` need to be modified
 ****************/
	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

	// fprintf(stderr, "&&&& sink id: %d;  graphno: %d \n", sink->nodeid, graphno);
/****************
 **  `no2gnode` need to be modified
 ****************/

	if(mergeMode) { // I might have a bunch of sink's parents that are not linked to sink
		for(int i=0;i<sink->parent.Count();i++) {
			CGraphnode *node=no2gnode[s][g][sink->parent[i]];
			node->child.Add(sink->nodeid);
		}
	}

	// fprintf(stderr,"This graph has %d nodes and %d edges and starts at lastpos=%d\n",graphno,edgeno,graphno);
	lastgpos=graphno; // nodes are from 0 to graphno-1, so the first "available" position in GBitVec is graphno

	/*****************************
	 ** Step 6: 'CTransfrag' constructor
	 ** 	now I can create the future transfrags because I know graphno
	 *****************************/
	for(int i=0;i<futuretr.Count();i+=3) {
		// add links between node and sink
		int n1=int(futuretr[i]);
		if(n1 >=0 ) {
			int n2=int(futuretr[i+1]);
			GBitVec trpat(graphno+edgeno);
			trpat[n1]=1;
			GVec<int> nodes;
			if(n2<0) {
				CGraphnode *node=no2gnode[s][g][n1];

				int n3=n1+1;
				if(n3<sink->nodeid) {
					CGraphnode *prevnode=node;
					CGraphnode *nextnode=no2gnode[s][g][n3]; // this should be okay
					uint dist=nextnode->len()-1;
					bool keeptr=true;
					while(nextnode->start==prevnode->end+1 && dist<longintronanchor) {
						for(int c=0;c<nextnode->child.Count();c++)
							if(nextnode->child[c]==sink->nodeid) {
								// do not continue with the futuretr
								keeptr=false;
								break;
							}
						n3++;
						if(n3==sink->nodeid) {
							// no link to sink found -> can continue with futuretr
							break;
						}
						prevnode=nextnode;
						nextnode=no2gnode[s][g][n3];
						dist+=nextnode->len();
					}
					if(!keeptr) continue; // found link to sink; no need to keep futuretr
				}

				node->child.Add(sink->nodeid);
				// add node to sink transfrag
				trpat[graphno-1]=1;

				int key=edge(n1,graphno-1,graphno);
				int *pos=gpos[s][g][key];
				if(pos!=NULL) trpat[*pos]=1;
				else {
					gpos[s][g].Add(key,lastgpos);
					trpat[lastgpos]=1;
					lastgpos++;
				}

				nodes.Add(n1);
				nodes.Add(sink->nodeid);
			}
			else {
				trpat[n2]=1;
				int key=edge(n1,n2,graphno);
				int *pos=gpos[s][g][key];
				if(pos!=NULL) trpat[*pos]=1;
				else {
					gpos[s][g].Add(key,lastgpos);
					trpat[lastgpos]=1;
					lastgpos++;
				}

				nodes.cAdd(n1);
				nodes.Add(n2);
			}
			CTransfrag *tr=new CTransfrag(nodes,trpat,futuretr[i+2]);

			// /*
			{ // DEBUG ONLY
				fprintf(stderr,"Add future transfrag[%d][%d]= %d with %d nodes n1=%d n2=%d graphno=%d, abundance=%f and pattern",s,g,transfrag[s][g].Count(),tr->nodes.Count(),n1,n2,graphno,futuretr[i+2]);
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			// */

			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
		}
	}

	/*****************************
	 ** Step 7: 'traverse_dfs_multi' function
	 ** 	finished reading bundle -> now create the parents' and children's patterns
	 *****************************/
	GVec<bool> visit;
	visit.Resize(graphno);
	GBitVec parents(graphno+edgeno);

	fprintf(stderr,"traverse graph[%d][%d] now with %d nodes, %d edges and lastgpos=%d....\n",s,g,graphno,edgeno,lastgpos);//edgeno=0;
	traverse_dfs_multi(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	fprintf(stderr,"done traversing with edgeno=%d lastgpos=%d\n",edgeno,lastgpos);

	// /*
	{ //DEBUG ONLY
		fprintf(stderr,"after traverse:\n");
		for(int i=0;i<graphno;i++) {
			fprintf(stderr,"Node %d with parents:",i);
			for(int p=0;p<no2gnode[s][g][i]->parent.Count();p++) fprintf(stderr," %d",no2gnode[s][g][i]->parent[p]);
			fprintf(stderr," and children:");
			for(int c=0;c<no2gnode[s][g][i]->child.Count();c++) fprintf(stderr," %d",no2gnode[s][g][i]->child[c]);
			fprintf(stderr,"\n");
		}
	}
	// */

	// delete variables created here, like e.g. ends; do I need to delete the GVec<int> elements created too?
	ends.Clear();

	return(graphno);
}






CTreePat *construct_treepat_multi(int gno, GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag) {

	/**********************
	 ** create root CTreePat first
	 **********************/
	CTreePat *root=new CTreePat(0,gno-1); // if links from source to nodes are desired source==1 and all nodes are treated as +1

	/**********************
	 ** now construct all child CTreePat's
	 **********************/
	// fprintf(stderr,"There are %d transfrags\n",transfrag.Count());
	for(int t=0;t<transfrag.Count();t++)
		if(transfrag[t]->nodes[0]){ // don't include transfrags from source -> not needed
			CTreePat *tree=root;
			int m=0; // previous node in pattern that was set in pattern
			for(int n=1;n<gno;n++){
				if(transfrag[t]->pattern[n]) {
					CTreePat *child;
					if(m) { // there is a node m that was seen before
						int *pos=gpos[edge(m,n,gno)];
						if(pos && transfrag[t]->pattern[*pos]) // there is an edge between m and n
							child=tree->settree(gno-1-m+n-m-1,n,2*(gno-n-1));
						else child=tree->settree(n-m-1,n,2*(gno-n-1));
					}
					else // this is the root tree
						child=tree->settree(n-1,n,2*(gno-n-1));
					tree=child;
					m=n;
				}
			}
			tree->tr=transfrag[t];
		}

	return(root);
}









int build_graphs_multi(BundleData* bdata, UniSpliceGraphGp* uni_splice_graphGp) {
	int refstart = bdata->start;
	int refend = bdata->end+1;
	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GPVec<GffObj>& guides = bdata->keepguides;
	GVec<float>* bpcov = bdata->bpcov; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles
	GList<CPrediction>& pred = bdata->pred;
	// form groups on strands: all groups below are like this: 0 = negative strand; 1 = unknown strand; 2 = positive strand
	GPVec<CGroup> group;
	CGroup *currgroup[3]={NULL,NULL,NULL}; // current group of each type
	CGroup *startgroup[3]={NULL,NULL,NULL}; // start group of each type
	int color=0; // next color to assign
	GVec<int> merge; // remembers merged groups
	GVec<int> equalcolor; // remembers colors for the same bundle
	GVec<int> *readgroup=new GVec<int>[readlist.Count()]; // remebers groups for each read; don't forget to delete it when no longer needed
	// parameter -e ; for mergeMode includes estimated coverage sum in the merged transcripts
	GVec<int> guidepred; // for eonly keeps the prediction number associated with a guide
	GArray<GEdge> guideedge; // 0: negative starts; 1 positive starts
	/*
	GPVec<GPtFeature>& feature = bdata->ptfs; // these are point features (confirmed starts/stops)

	for(int i=0;i<feature.Count();i++) {
		if(feature[i]->ftype==GPFT_TSS)
			fprintf(stderr,"TSS at position %d on strand %d\n",feature[i]->coord,feature[i]->strand);
		if(feature[i]->ftype==GPFT_CPAS)
			fprintf(stderr,"CPAS at position %d on strand %d\n",feature[i]->coord,feature[i]->strand);
	}
	*/

	fprintf(stderr,"build_graphs with %d guides\n",guides.Count());

	/*****************************
	 ** Step 1: this part is for setting guides for introns are covered by at least one read. 
	 *****************************/
	if(guides.Count()) {

		guideedge.setSorted(true);
		guideedge.setUnique(true);

		//if(eonly)
		for(int g=0;g<guides.Count();g++) {
			guidepred.cAdd(-1);
			bool covered=true;
			RC_TData* tdata=(RC_TData*)(guides[g]->uptr);
			if(longreads || mixedMode) {
				for(int i=1;i<guides[g]->exons.Count();i++) {
					char s=0; // unknown strand
					if(guides[g]->strand=='+') s=1; // guide on positive strand
					else if(guides[g]->strand=='-') s=-1; // guide on negative strand
					CJunction jn(guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s);
					int oidx=-1;
					if (!junction.Found(&jn, oidx)) {
						covered=false;
						break;
					}
				}
			} else {
				for(int i=0;i<tdata->t_introns.Count();i++) {
					if(!tdata->t_introns[i]->rcount) {
						covered=false;
						break;
					}
				}
			}

			if(covered) {
				// 1 if used by read bundles (present in keepguides) 
				// 2 if all introns are covered by at least one read
				// 3 if it is stored to be printed
				tdata->in_bundle=2;
				int s=-1; // unknown strand
				if(guides[g]->strand=='+') s=1; // guide on positive strand
				else if(guides[g]->strand=='-') s=0; // guide on negative strand

				fprintf(stderr,"Look to add guide g=%d start %d-%d and end %d-%d on strand %d\n",g,guides[g]->start,guides[g]->exons[0]->end,guides[g]->end,guides[g]->exons.Last()->start,s);

				int uses=s;
				if(s<0) uses=0;
				// guide edge
				GEdge ge(guides[g]->start,guides[g]->exons[0]->end,uses);
				int idx=guideedge.IndexOf(ge);

				fprintf(stderr,"look for ge(%d,%d,%d) => start idx=%d\n",ge.val,ge.endval,ge.strand,idx);

				if(idx<0) guideedge.Add(ge);
				else if(guideedge[idx].endval>guides[g]->exons[0]->end) guideedge[idx].endval=guides[g]->exons[0]->end;
				if(s<0) {
					ge.strand=1;
					idx=guideedge.IndexOf(ge);
					if(idx<0) guideedge.Add(ge);
					else if(guideedge[idx].endval>guides[g]->exons[0]->end) guideedge[idx].endval=guides[g]->exons[0]->end;
				}

				ge.val=guides[g]->end;
				ge.endval=guides[g]->exons.Last()->start;
				idx=guideedge.IndexOf(ge);

				fprintf(stderr,"look for ge(%d,%d,%d) => end idx=%d\n",ge.val,ge.endval,ge.strand,idx);

				if(idx<0) guideedge.Add(ge);
				else if(guideedge[idx].endval<guides[g]->exons.Last()->start) guideedge[idx].endval=guides[g]->exons.Last()->start;
				if(s<0) {
					ge.strand=0; // ge.strand was set as 1 before
					idx=guideedge.IndexOf(ge);
					if(idx<0) guideedge.Add(ge);
					else if(guideedge[idx].endval<guides[g]->exons.Last()->start) guideedge[idx].endval=guides[g]->exons.Last()->start;
				}
			}
		}
	}

	/*****************************
	 ** Step 2: this part is for adjusting leftsupport and rightsupport when considering all junctions that start at a given point
	 ** 	sort junctions -> junctions are sorted already according with their start, but not their end
	 *****************************/
	GList<CJunction> ejunction(junction);
	ejunction.setFreeItem(false);
	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);

	uint start=0;
	uint end=0;
	double leftsupport[2]={0,0};    // should be strand based
	double rightsupport[2]={0,0};
	bool higherr=false;
	char leftcons=-1;
	char rightcons=-1;
	for(int i=0;i<junction.Count();i++) {

		//fprintf(stderr,"check junction:%d-%d:%d leftsupport=%f rightsupport=%f nm=%f nreads=%f\n",junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->leftsupport,junction[i]->rightsupport,junction[i]->nm,junction[i]->nreads);

		if((!higherr || mixedMode) && junction[i]->strand && junction[i]->nm==junction[i]->nreads && !junction[i]->guide_match) {
			higherr=true;
			if(mixedMode) {
				int j=i-1;
				while(j>=0 && junction[j]->start+sserror>junction[i]->start) {
					if(junction[j]->strand && junction[i]->strand!=junction[j]->strand && abs((int)(junction[j]->end-junction[i]->end))<(int)sserror && junction[j]->nm<junction[j]->nreads) {
						junction[i]->strand = 0;
						break;
					}
					j--;
				}
				if(junction[i]->strand) {
					j=i+1;
					while(j<junction.Count() && junction[j]->start-sserror<junction[i]->start) {
						if(junction[j]->strand && junction[i]->strand!=junction[j]->strand && abs((int)(junction[j]->end-junction[i]->end))<(int)sserror && junction[j]->nm<junction[j]->nreads) {
							junction[i]->strand = 0;
							break;
						}
						j++;
					}
				}
			}
		}

		if(junction[i]->start!=start) {  // new junction starting at i

			int j=i-1;
			while(j>=0 && junction[j]->start==start) {
				junction[j]->consleft=leftcons;
				if(!junction[j]->guide_match && !leftcons && junction[j]->nreads_good<DROP/ERROR_PERC) junction[j]->strand=0;
				if(junction[j]->strand) junction[j]->leftsupport=leftsupport[(1+junction[j]->strand)/2];
				j--;
			}


			j=i+1; // check if there is the same junction further ahead first
			if(junction[i]->strand) while(j<junction.Count() && junction[j]->start==junction[i]->start && junction[j]->end==junction[i]->end) {
				if(junction[j]->strand && junction[i]->strand!=junction[j]->strand) {
					//possible missaligned junction --> delete junction strand
					if(junction[i]->nreads>junction[j]->nreads && junction[i]->nreads_good>junction[j]->nreads_good) junction[j]->strand=0;
					if(junction[i]->nreads<junction[j]->nreads && junction[i]->nreads_good<junction[j]->nreads_good) junction[i]->strand=0;
				}
				j++;
			}

			leftsupport[0]=0;
			leftsupport[1]=0;
			if(junction[i]->strand) leftsupport[(1+junction[i]->strand)/2]=junction[i]->leftsupport; // I might have deleted the i junction above
			start=junction[i]->start;

			if(bdata->gseq) {
				if(junction[i]->strand>0) {
					if((bdata->gseq[junction[i]->start+1-refstart]=='g'||bdata->gseq[junction[i]->start+1-refstart]=='G') &&
							(bdata->gseq[junction[i]->start+2-refstart]=='T'||bdata->gseq[junction[i]->start+2-refstart]=='t'||
									bdata->gseq[junction[i]->start+2-refstart]=='C'||bdata->gseq[junction[i]->start+2-refstart]=='c')) {
						junction[i]->consleft=1;
						leftcons=1;
					}
					else {
						junction[i]->consleft=0;
						leftcons=0;
					}
				}
				else if(junction[i]->strand<0) {
					if((bdata->gseq[junction[i]->start+1-refstart]=='c'||bdata->gseq[junction[i]->start+1-refstart]=='C') &&
							(bdata->gseq[junction[i]->start+2-refstart]=='T'||bdata->gseq[junction[i]->start+2-refstart]=='t')) {
						junction[i]->consleft=1;
						leftcons=1;
					}
					else {
						junction[i]->consleft=0;
						leftcons=0;
					}
				}
				//fprintf(stderr,"junction:%d-%d:%d is %c%c-%c%c\n",junction[i]->start,junction[i]->end,junction[i]->strand,bdata->gseq[junction[i]->start+1-refstart],bdata->gseq[junction[i]->start+2-refstart],bdata->gseq[junction[i]->end-2-refstart],bdata->gseq[junction[i]->end-1-refstart]);
			}
		}
		else if(junction[i]->strand) leftsupport[(1+junction[i]->strand)/2]+=junction[i]->leftsupport;
		//fprintf(stderr,"leftsupport[%d]=%f\n",(1+junction[i]->strand)/2,leftsupport[(1+junction[i]->strand)/2]);


		//fprintf(stderr,"check ejunction:%d-%d:%d leftsupport=%f rightsupport=%f nm=%f nreads=%f\n",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->leftsupport,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);
		if(ejunction[i]->end!=end) { // I do not check if I deleted the junction here for support

			int j=i-1;
			while(j>=0 && ejunction[j]->end==end) {
				ejunction[j]->consright=rightcons;
				if(!ejunction[j]->guide_match && !rightcons && ejunction[j]->nreads_good<DROP/ERROR_PERC) ejunction[j]->strand=0;
				if(ejunction[j]->strand) ejunction[j]->rightsupport=rightsupport[(1+ejunction[j]->strand)/2];
				j--;
			}
			// I do not check here for possible missalignments
			rightsupport[0]=0;
			rightsupport[1]=0;
			if(ejunction[i]->strand) rightsupport[(1+ejunction[i]->strand)/2]=ejunction[i]->rightsupport;
			end=ejunction[i]->end;

			if(bdata->gseq) {
				if(ejunction[i]->strand>0) {
					if((bdata->gseq[ejunction[i]->end-2-refstart]=='a'||bdata->gseq[ejunction[i]->end-2-refstart]=='A') &&
							(bdata->gseq[ejunction[i]->end-1-refstart]=='G'||bdata->gseq[ejunction[i]->end-1-refstart]=='g')) {
						ejunction[i]->consright=1;
						rightcons=1;
					}
					else {
						ejunction[i]->consright=0;
						rightcons=0;
					}
				}
				else if(ejunction[i]->strand<0) {
					if((bdata->gseq[ejunction[i]->end-1-refstart]=='C'||bdata->gseq[ejunction[i]->end-1-refstart]=='c') &&
							(bdata->gseq[ejunction[i]->end-2-refstart]=='A'||bdata->gseq[ejunction[i]->end-2-refstart]=='a'||
									bdata->gseq[ejunction[i]->end-2-refstart]=='G'||bdata->gseq[ejunction[i]->end-2-refstart]=='g')) {
						ejunction[i]->consright=1;
						rightcons=1;
					}
					else {
						ejunction[i]->consright=0;
						rightcons=0;
					}
				}
				//fprintf(stderr,"ejunction:%d-%d:%d is %c%c-%c%c\n",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,bdata->gseq[ejunction[i]->start+1-refstart],bdata->gseq[ejunction[i]->start+2-refstart],bdata->gseq[ejunction[i]->end-2-refstart],bdata->gseq[ejunction[i]->end-1-refstart]);
			}

		}
		else if(ejunction[i]->strand) rightsupport[(1+ejunction[i]->strand)/2]+=ejunction[i]->rightsupport;
		//fprintf(stderr,"rightsupport[%d]=%f\n",(1+ejunction[i]->strand)/2,rightsupport[(1+ejunction[i]->strand)/2]);
	}
	// end adjusting leftsupport and rightsupport

	//fprintf(stderr,"junction support computed\n");

	/*****************************
	 ** Step 3: there are some reads that contain very bad junctions -> need to find better closest junctions
	 *****************************/
	if(higherr) { 
		uint juncsupport=junctionsupport;
		if(longreads)
			juncsupport=sserror;
		else if(mixedMode) juncsupport=sserror/DROP;
		//fprintf(stderr,"In higherr!\n");
		GVec<int> jstarts; // good starting junctions
		GVec<int> jends; // good ending junctions
		if(viral) {
			for(int i=1;i<junction.Count();i++) { // junction is sorted based on start

				if(junction[i]->strand && junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					if(junction[i]->nreads_good>=0 && (junction[i]->nreads_good<1.25*junctionthr || !good_junc(*junction[i],refstart,bpcov))) { // threshold for bad junctions is higher; (should I also add that too short junctions not to be accepted?)
						//junction[i]->strand=0; // just delete junction if it's low count
						junction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}

					int j=i-1;
					while(j>0 && junction[i]->start-junction[j]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand && junction[j]->nreads_good>=0 && junction[i]->nreads<junction[j]->nreads &&
								abs((int)junction[i]->end-(int)junction[j]->end)<(int)juncsupport) { // j was not elminated and is better
							if(junction[i]->nreads_good<0) { // i was eliminated before
								int k=-junction[i]->nreads_good;
								if(junction[k]->nreads<junction[j]->nreads) junction[i]->nreads_good=-j;
							}
							else {
								junction[i]->nreads_good=-j;
							}
						}
						j--;
					}
				}
			}
			for(int i=junction.Count()-1;i>0;i--) {
				if(junction[i]->strand && junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					int j=i+1;
					while(j<junction.Count() && junction[j]->start-junction[i]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand && junction[j]->nreads_good>=0 && junction[i]->nreads<junction[j]->nreads &&
								abs((int)junction[i]->end-(int)junction[j]->end)<(int)juncsupport) { // j was not elminated and is better
							if(junction[i]->nreads_good<0) { // i was eliminated before
								int k=-junction[i]->nreads_good;
								if(junction[k]->nreads<junction[j]->nreads) junction[i]->nreads_good=-j;
							}
							else {
								junction[i]->nreads_good=-j;
							}
						}
						j++;
					}
				}
			}
		} else {

			float tolerance=1-ERROR_PERC;

			// strand based version
			for(int i=1;i<junction.Count();i++) { // junction is sorted based on start

				//fprintf(stderr,"junct[%d]:%d-%d:%d lefttsupport=%f nm=%f mm=%f nreads=%f nreads_good=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->leftsupport,junction[i]->nm,junction[i]->mm,junction[i]->nreads,junction[i]->nreads_good);
				if(junction[i]->strand) {
					if(junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					if(junction[i]->nreads_good>=0 && junction[i]->nreads_good<1.25*junctionthr) { // threshold for bad junctions is higher; (should I also add that too short junctions not to be accepted?)
						//junction[i]->strand=0; // just delete junction if it's low count
						junction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}

					int j=i-1;
					float support=0;
					bool searchjunc=true;
					bool reliable=false;
					//if(j>=0) fprintf(stderr,"...start at junct:%d-%d:%d leftsupport=%f dist=%d\n",junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport,junction[i]->start-junction[j]->start);
					while(j>0 && junction[i]->start-junction[j]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand) {
							if(junction[j]->start==junction[i]->start) { // found a junction with the same start -> I have already searched it if it's bad
								if(junction[j]->nreads<0) {
									junction[i]->nreads=junction[j]->nreads;
									searchjunc=false;
								}
								break;
							}
							else if(junction[j]->guide_match || junction[j]->nm<junction[j]->nreads) { // nearby junction is much more reliable
								if(junction[j]->leftsupport>junction[i]->leftsupport*tolerance) { // the good junction is close enough
									reliable=true;
									junction[i]->nreads=-j;
									support=junction[j]->leftsupport;
									break;
								}
							}
							else if(junction[j]->leftsupport>support && junction[i]->start-junction[j]->start<sserror && junction[j]->leftsupport*tolerance>junction[i]->leftsupport) {
								//fprintf(stderr,"...1 compare to [%d]:%d-%d:%d leftsupport=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport);
								junction[i]->nreads=-j;
								support=junction[j]->leftsupport;
							}
						}
						j--;
					}
					if(searchjunc) {
						j=i+1;
						int dist=juncsupport;
						if(junction[i]->nreads<0) {
							dist=junction[i]->start-junction[abs((int)junction[i]->nreads)]->start;
						}
						//if(j<junction.Count()) fprintf(stderr,"...start at junct:%d-%d:%d leftsupport=%f dist=%d\n",junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport,junction[j]->start-junction[i]->start);
						while(j<junction.Count() && junction[j]->start-junction[i]->start<juncsupport) {
							if(junction[j]->strand==junction[i]->strand && junction[j]->start!=junction[i]->start) {
								int d=(int)(junction[j]->start-junction[i]->start);
								if(junction[j]->guide_match || junction[j]->nm<junction[j]->nreads) {
									if((d<dist || (d==dist && junction[j]->leftsupport>support)) && junction[j]->leftsupport>junction[i]->leftsupport*tolerance) {
										junction[i]->nreads=-j;
										support=junction[j]->leftsupport;
										break;
									}
								}
								else if(!reliable && junction[j]->leftsupport>support && (uint)d<sserror && junction[j]->leftsupport*tolerance>junction[i]->leftsupport) { // junction is not best within window
									//fprintf(stderr,"...2 compare to [%d]:%d-%d:%d leftsupport=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport);
									junction[i]->nreads=-j;
									support=junction[j]->leftsupport;
								}
							}
							j++;
						}
					}
				}
					else if(mixedMode && junction[i]->nm<junction[i]->nreads) {
						jstarts.Add(i);
					}
				}
				//fprintf(stderr,"ejunct[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);

				if(ejunction[i]->strand) {
					if(ejunction[i]->nm && !ejunction[i]->guide_match && ejunction[i]->nm>=ejunction[i]->nreads) { // this is a bad junction -> check if it's maximal
					if(ejunction[i]->nreads_good>=0 && ejunction[i]->nreads_good<1.25*junctionthr) { // threshold for bad junctions is higher
						//ejunction[i]->strand=0;
						ejunction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}
					int j=i-1;
					float support=0;
					bool searchjunc=true;
					bool reliable=false;
					//if(j>=0) fprintf(stderr,"...start at junct:%d-%d:%d rightsupport=%f dist=%d\n",ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport,ejunction[i]->end-ejunction[j]->end);
					while(j>0 && ejunction[i]->end-ejunction[j]->end<juncsupport) {
						if(ejunction[j]->strand==ejunction[i]->strand) {
							if(ejunction[j]->end==ejunction[i]->end) {
								if(ejunction[j]->nreads_good<0) {
									ejunction[i]->nreads_good=ejunction[j]->nreads_good;
									searchjunc=false;
								}
								break;
							}
							else if(ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads) { // nearby junction is much more reliable
								if(ejunction[j]->rightsupport>ejunction[i]->rightsupport*tolerance) { // the good junction is close enough
									reliable=true;
									ejunction[i]->nreads_good=-j;
									support=ejunction[j]->rightsupport;
									break;
								}
							}
							else if(ejunction[j]->rightsupport>support && ejunction[i]->end-ejunction[j]->end < sserror && ejunction[j]->rightsupport*tolerance>ejunction[i]->rightsupport) {
								//fprintf(stderr,"...1 compare to [%d]:%d-%d:%d rightsupport=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport);
								ejunction[i]->nreads_good=-j;
								support=ejunction[j]->rightsupport;
							}
						}
						j--;
					}
					if(searchjunc) {
						j=i+1;
						int dist=juncsupport;
						if(ejunction[i]->nreads_good<0) {
							dist=ejunction[i]->end-ejunction[abs((int)ejunction[i]->nreads_good)]->end;
						}
						//if(j<junction.Count()) fprintf(stderr,"...start at junct:%d-%d:%d rightsupport=%f dist=%d\n",ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport,ejunction[j]->end-ejunction[i]->end);
						while(j<junction.Count() && ejunction[j]->end-ejunction[i]->end<juncsupport) {
							if(ejunction[j]->strand==ejunction[i]->strand && ejunction[j]->end!=ejunction[i]->end) {
								int d=ejunction[j]->end-ejunction[i]->end;
								if(ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads) {
									if((d<dist || (d==dist && ejunction[j]->rightsupport>support)) && ejunction[j]->rightsupport>ejunction[i]->rightsupport*tolerance) {
										ejunction[i]->nreads_good=-j;
										support=ejunction[j]->rightsupport;
										break;
									}
								}
								else if((!reliable && ejunction[j]->rightsupport>support && ejunction[j]->end-ejunction[i]->end < sserror && ejunction[j]->rightsupport*tolerance>ejunction[i]->rightsupport) || ((int)(ejunction[j]->end-ejunction[i]->end)<dist && (ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads))) {
									//fprintf(stderr,"...2 compare to [%d]:%d-%d:%d rightsupport=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport);
									ejunction[i]->nreads_good=-j;
									support=ejunction[j]->rightsupport;
								}
							}
							j++;
						}
					}
				}
					else if(mixedMode && ejunction[i]->nm<ejunction[i]->nreads) {
						jends.Add(i);
					}
				}
			}
		}
		if(mixedMode) { // check if there are junctions inside bigger junctions that can form small exons
			int si=0;
			for(int ei=0;ei<jends.Count();ei++) {
				while(si<jstarts.Count() && junction[jstarts[si]]->start<=ejunction[jends[ei]]->end) si++;
				int k=si;
				while(k<jstarts.Count() && junction[jstarts[k]]->start-ejunction[jends[ei]]->end<SMALL_EXON) {
					if(junction[jstarts[k]]->strand == ejunction[jends[ei]]->strand) {
						CJunction jn(ejunction[jends[ei]]->start,junction[jstarts[k]]->end,junction[jstarts[k]]->strand);
						int oidx=-1;
						if (junction.Found(&jn, oidx) && junction[oidx]->nm>=junction[oidx]->nreads) { // candidate junction for deletion
							junction[oidx]->strand=0;
							break;
						}
					}
					k++;
				}
			}
		}
	} //if higherr

	/*
	{ // DEBUG ONLY
		for(int i=0;i<junction.Count();i++) {
			if(junction[i]->guide_match) fprintf(stderr,"G");
			if(junction[i]->strand && junction[i]->nreads>0 && junction[i]->nreads_good>0 && junction[i]->mm>=0) fprintf(stderr,"***");
			//if(junction[i]->strand && junction[i]->mm>=0)
				fprintf(stderr,"Junction[%d] %d-%d:%d has nm=%f mm=%f nreads=%f nreads_good=%f leftsupport=%f and rightsupport=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,
					junction[i]->nm,junction[i]->mm,junction[i]->nreads,
					junction[i]->nreads_good,junction[i]->leftsupport,junction[i]->rightsupport);
		}
	}
	*/
	//int **readgroup = new int*[readlist.Count()];
/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	/*****************************
	 ** Step 4: junctions filtering & sort
	 *****************************/
	//float fraglen=0;
	//uint fragno=0;
	//GHash<bool> boundaryleft;
	//GHash<bool> boundaryright;
	GIntHash<bool> boundaryleft;
	GIntHash<bool> boundaryright;

	bool resort=false;
	int njunc=junction.Count();
	for (int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		// version that does not adjust read but discards it completely
		bool keep=true;
		int i=0;

		/*fprintf(stderr,"check read[%d]:%d-%d:%d refstart=%d w/exons:",n,readlist[n]->start,readlist[n]->end,readlist[n]->strand,refstart);
		for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
		fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);
		i=0;*/


		while(i<rd.juncs.Count()) {
			CJunction& jd=*(rd.juncs[i]);
			//fprintf(stderr, " read junc %d-%d:%d nreads=%f nreads_good=%f nm=%f support=%f,%f between exons:%d-%d and %d-%d\n", jd.start, jd.end, jd.strand,jd.nreads,jd.nreads_good,jd.nm,jd.leftsupport,jd.rightsupport,rd.segs[i].start,rd.segs[i].end,rd.segs[i+1].start,rd.segs[i+1].end);

			bool changeright=jd.nreads_good<0;
			bool changeleft=jd.nreads<0;
			if(viral && changeright) {
				changeleft=true;
				jd.nreads=jd.nreads_good;
			}

			if(jd.strand) {
				if(!good_junc(jd,refstart,bpcov)) jd.mm=-1; // bad junction
				if(changeleft || changeright) jd.strand=0;
				if(jd.mm<0) jd.strand=0; // I can't believe the junction because it's too low
			}

			if(!jd.strand) { // found a bad junction
				//instead of deleting read -> just delete junction and make sure that add_read_to_group adds read to multiple groups
				//jd.nreads-=rd.read_count;

				//if(jd.mm>=0 && (changeleft || changeright)) { // not sure why jd.mm should be positive here since I am searching the junction against the others anyway
				if(changeleft || changeright) {
					uint newstart=rd.segs[i].end;
					uint newend=rd.segs[i+1].start;
					int jk=-1;
					int ek=-1;
					if(changeleft) {
						jk=abs(int(jd.nreads));
						if(junction[jk]->nreads<0) { // not a valid left support -> delete junction
							newstart=rd.segs[i].start-1;
						}
						else {
							newstart=junction[jk]->start; // junction jk is good
							//fprintf(stderr,"junction has newstart=%d from junction[jk=%d]\n",newstart,jk);
							if(junction[jk]->strand) jk=-1;
						}
					}
					if(changeright) {
						ek=abs(int(jd.nreads_good));
						if(viral) {
							if(junction[ek]->nreads_good<0) {
								newend=rd.segs[i+1].end+1;
							}
							else {
								newend=junction[ek]->end;
								//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
							}
						}
						else {
							if(ejunction[ek]->nreads_good<0) {
								newend=rd.segs[i+1].end+1;
							}
							else {
								newend=ejunction[ek]->end;
								//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
								if(ejunction[ek]->strand) ek=-1;
							}
						}
					}
					//if(jd.mm>=0 && newstart>=rd.segs[i].start && newend<=rd.segs[i+1].end) { // junction inside read boundaries
					if(newstart>=rd.segs[i].start && newend<=rd.segs[i+1].end && newstart<=newend) { // junction inside read boundaries
						bool searchjunc=true;
						bool addjunction=true;
						rd.segs[i].end=newstart;   // adjust start
						rd.segs[i+1].start=newend; // adjust end
						if(jd.mm>=0) {
							if(jk>0) { // search junctions
								int k=jk;
								while(k>0 && junction[k]->start==newstart) {
								   if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
								   if(rd.strand==junction[k]->strand && junction[k]->end==newend) {
									rd.juncs.Put(i,junction[k]);
										searchjunc=false;
										break;
									}
									k--;
								}
								if(searchjunc) {
									k=jk+1;
									while(k<junction.Count() && junction[k]->start==newstart) {
										if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
										if(rd.strand==junction[k]->strand && junction[k]->end==newend) {
											rd.juncs.Put(i,junction[k]);
											searchjunc=false;
											break;
										}
										k++;
									}
								}
								if(searchjunc) { // I did not find junction -> I need to create a new one
									// first check if I already added such a junction
									for(k=njunc;k<junction.Count();k++) {
										if(rd.strand==junction[k]->strand && junction[k]->start==newstart && junction[k]->end==newend) {
											rd.juncs.Put(i,junction[k]);
											searchjunc=false;
											break;
										}
									}
									if(searchjunc && addjunction) {
										if(!resort) {
											junction.setSorted(false);
											ejunction.setSorted(false);
											resort=true;
										}
										//fprintf(stderr,"Add new junction:%d-%d:%d at position %d njunc=%d\n",newstart,newend,rd.strand,junction.Count(),njunc);
										CJunction *junc=new CJunction(newstart,newend,rd.strand);
										junction.Add(junc);
										ejunction.Add(junc);
										rd.juncs.Put(i,junc);
									}
								}
							}
							else if(ek>0) { // ek>0 => search ejunctions; ek>0 because I have either changeleft or changeright
								int k=ek;
								while(k>0 && ejunction[k]->end==newend) {
									if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
									if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart) {
										rd.juncs.Put(i,ejunction[k]);
										searchjunc=false;
										break;
									}
									k--;
								}
								if(searchjunc) {
									k=ek+1;
									while(k<ejunction.Count() && ejunction[k]->end==newend) {
										if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
										if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart) {
											rd.juncs.Put(i,ejunction[k]);
											searchjunc=false;
											break;
										}
										k++;
									}
								}
								if(searchjunc && addjunction) { // I did not find junction -> I need to create a new one
									// first check if I already added such a junction
									for(k=njunc;k<ejunction.Count();k++) {
										if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart && ejunction[k]->end==newend) {
											rd.juncs.Put(i,ejunction[k]);
											searchjunc=false;
											break;
										}
									}
									if(searchjunc) {
										if(!resort) {
											junction.setSorted(false);
											ejunction.setSorted(false);
											resort=true;
										}
										//fprintf(stderr,"Add new junction:%d-%d:%d at position %d njunc=%d\n",newstart,newend,rd.strand,junction.Count(),njunc);
										CJunction *junc=new CJunction(newstart,newend,rd.strand);
										junction.Add(junc);
										ejunction.Add(junc);
										rd.juncs.Put(i,junc);
									}
								}
							}
						}
					}
					//fprintf(stderr, "read[%d] adjusted to junction:%d-%d\n",n,rd.segs[i].end,rd.segs[i+1].start);
				}

				// because read might be poorly mapped I also have to unpair it
				if(keep) { // this is the first time I unpair read -> remove strand if pair has single exon?
					for(int p=0;p<rd.pair_idx.Count();p++) {
						int np=rd.pair_idx[p];
						if(np>-1) {
							rd.pair_idx[p]=-1;
							for(int j=0;j<readlist[np]->pair_idx.Count();j++) // also unpair read np to n
								if(readlist[np]->pair_idx[j]==n) {
									readlist[np]->pair_idx[j]=-1;
									break;
								}
							// remove strand for pair read if single exon ? is strand assigned before here ?
							if(!readlist[np]->juncs.Count()) readlist[np]->strand=0;
						}
					}
					keep=false;
				}
			}
			else if(guides.Count()){ // need to remember boundary
				bool exist=true;
		    	//GStr bs((int)jd.start);
		    	if(!boundaryleft[jd.start]) boundaryleft.Add(jd.start,exist);
		    	//GStr be((int)jd.end);
		    	if(!boundaryright[jd.end]) boundaryright.Add(jd.end,exist);
			}
			i++;

		}

		if(!keep) { // read has bad junctions -> check if I need to unstrand it
			if(rd.strand) {
				bool keepstrand=false;
				for(int j=0;j<rd.juncs.Count();j++) if(rd.juncs[j]->strand) { keepstrand=true; break; }
				if(!keepstrand) rd.strand=0;
			}
		}
		else if(!rd.juncs.Count() && rd.strand) { // check if I need to unstrand current read due to future mapping
			if(rd.pair_idx.Count()) {
				bool keepstrand=false;
				for(int p=0;p<rd.pair_idx.Count();p++) {
					int np=rd.pair_idx[p];
					if(np>-1 && n<np) {
						if(!readlist[np]->juncs.Count()) {
							if(readlist[np]->strand==rd.strand) {
								keepstrand=true;
								break;
							}
						}
						else {
							for(int j=0;j<readlist[np]->juncs.Count();j++) {
								CJunction& jd=*(readlist[np]->juncs[j]);
								if(jd.strand && good_junc(jd,refstart,bpcov)) {
									keepstrand=true;
									break;
								}
							}
							if(!keepstrand) readlist[np]->strand=0;
						}
					}
					if(keepstrand) break;
				}
				if(!keepstrand) rd.strand=0;
			}
		}


		//if(rd.juncs.Count()) fprintf(stderr,"read[%d] keep=%d\n",n,keep);
		//if(rd.strand) fprintf(stderr,"read[%d] has strand %d\n",n,rd.strand);


		//if(keep) { // if it's a good read that needs to be kept


			/*fprintf(stderr,"add read %d:%d-%d w/count=%g for color=%d with npairs=%d\n",n,readlist[n]->start,readlist[n]->end,readlist[n]->read_count,color,readlist[n]->pair_idx.Count());
			fprintf(stderr,"add read[%d]:%d-%d:%d w/count=%g w/exons:",n,readlist[n]->start,readlist[n]->end,readlist[n]->strand,readlist[n]->read_count);
			for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
			fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);*/

			color=add_read_to_group(n,readlist,color,group,currgroup,startgroup,readgroup,equalcolor,merge);

			// count fragments
			if(!rd.unitig)
				bdata->frag_len+=rd.len*rd.read_count; // TODO: adjust this to work with FPKM for super-reads and Pacbio
			double single_count=rd.read_count;
			if(keep) for(int i=0;i<rd.pair_idx.Count();i++) {
				// I am not counting the fragment if I saw the pair before and it wasn't deleted
				if(rd.pair_idx[i]!=-1 && n>rd.pair_idx[i] && readlist[rd.pair_idx[i]]->nh) {// only if read is paired and it comes first in the pair I count the fragments
					single_count-=rd.pair_count[i];
				}
			}
			if(!rd.unitig && single_count>epsilon) {
				bdata->num_fragments+=single_count; // TODO: FPKM will not work for super-reads here because I have multiple fragments in
												    // a super-read -> I might want to re-estimate this from coverage and have some input for read length; or I might only use TPM
			}

			//fprintf(stderr,"now color=%d\n",color);
		//}
		//else { fprintf(stderr,"read[%d] is not kept\n",n);}
		//else clean_read_junctions(readlist[n]);
	}

	if(resort) {
		junction.setSorted(true);
		ejunction.setSorted(juncCmpEnd);
	}
	//fprintf(stderr,"fragno=%d fraglen=%g\n",fragno,fraglen);
	//if(fragno) fraglen/=fragno;

	/*****************************
	 ** Step 5: 'merge_fwd_groups' function
	 ** 	merge groups that are close together or __groups that are within the same exon of a reference gene__
	 *****************************/
	if(bundledist || (guides.Count() && !longreads)) {
		for(int sno=0;sno<3;sno++) {
			CGroup *lastgroup=NULL;
			CGroup *procgroup=startgroup[sno];
			while(procgroup!=NULL) {

				if(lastgroup) {

					//fprintf(stderr,"sno=%d lastgroup->end=%d procgroup->start=%d procgroup->end=%d\n",sno,lastgroup->end,procgroup->start,procgroup->end);

					//GStr bstart((int)lastgroup->end);
					//GStr bend((int)procgroup->start);
					if(!boundaryleft[lastgroup->end] && !boundaryright[procgroup->start] && (procgroup->start-lastgroup->end<=bundledist ||
			    				(guides.Count()  && guide_exon_overlap(guides,sno,lastgroup->end,procgroup->start)))) {

			    			//fprintf(stderr,"sno=%d merge groups btw %d and %d dist=%d\n",sno,lastgroup->end,procgroup->start,procgroup->start-lastgroup->end);

						merge_fwd_groups(group,lastgroup,procgroup,merge,equalcolor);
						procgroup=lastgroup->next_gr;
						continue;
					}
				}
				lastgroup=procgroup;
				procgroup=procgroup->next_gr;
			}
		}
	}

	if(guides.Count()) {
		boundaryleft.Clear();
		boundaryright.Clear();
	}
	/*
	{ // DEBUG ONLY
		fprintf(stderr,"%d groups created!\n",group.Count());
	    for(int sno=0;sno<3;sno++) {
	    	fprintf(stderr, "Groups on strand %d:\n",sno);
	    	CGroup *procgroup=startgroup[sno];
			while(procgroup!=NULL) {
	    		fprintf(stderr, " gr %d(%d,%.6f): %d-%d",procgroup->grid,procgroup->color,procgroup->cov_sum,procgroup->start,procgroup->end);
	    		procgroup=procgroup->next_gr;
	    	}
	    	fprintf(stderr,"\n");
	    }
	}
	*/
/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(after groups created):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	/*****************************
	 ** Step 6: form bundles here
     ** 	first do color assignment
	 *****************************/
	for (int i=0;i<3;i++) currgroup[i]=startgroup[i];
	CGroup *prevgroup[3]={NULL,NULL,NULL};

	GVec<int> eqposcol(true);
	GVec<int> eqnegcol(true);

	eqposcol.Resize(equalcolor.Count(),-1);
	eqnegcol.Resize(equalcolor.Count(),-1);


	// each unstranded group needs to remember what proportion of stranded group it overlaps so that it can distribute reads later on -> maybe I can do this in the following while?
	/*****************************
	 ** Step 7: 'set_strandcol' function
     ** 	set the color of strands
	 *****************************/
	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup); // gets the index of currgroup with the left most begining

		int grcol = currgroup[nextgr]->color;    // set smallest color for currgroup[$nextgr]

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;

		//fprintf(stderr,"group %d id=%d: %u-%u col=%d from col=%d\n",nextgr,currgroup[nextgr]->grid,currgroup[nextgr]->start,currgroup[nextgr]->end,grcol,prevcol);

		if(nextgr == 1) { // unknown strand group

			if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end+bundledist) { // overlaps previous negative group ; this needs bundledist
				//fprintf(stderr,"\tovlp to neg group: %u-%u\n",prevgroup[0]->start,prevgroup[0]->end);
				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > prevgroup[0]->start ? currgroup[nextgr]->start : prevgroup[0]->start;
				uint minend = currgroup[nextgr]->end < prevgroup[0]->end ? currgroup[nextgr]->end : prevgroup[0]->end;
				if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
				currgroup[nextgr]->neg_prop+=prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len();
			}

			while(currgroup[0]!=NULL && currgroup[nextgr]->start <= currgroup[0]->end+bundledist && currgroup[0]->start <= currgroup[nextgr]->end +bundledist) { // overlaps current negative strand group
				//fprintf(stderr,"\tovlp to neg group: %u-%u\n",currgroup[0]->start,currgroup[0]->end);

				int grcol = currgroup[0]->color;    // set smallest color for currgroup[$nextgr]
				while(equalcolor[grcol]!=grcol) {
					grcol=equalcolor[grcol];
				}
				currgroup[0]->color=grcol;
				currgroup[0]->neg_prop=1;

				set_strandcol(currgroup[nextgr],currgroup[0],currgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > currgroup[0]->start ? currgroup[nextgr]->start : currgroup[0]->start;
				uint minend = currgroup[nextgr]->end < currgroup[0]->end ? currgroup[nextgr]->end : currgroup[0]->end;
				if(minend<maxstart) minend=maxstart;
				currgroup[nextgr]->neg_prop+=currgroup[0]->cov_sum*(minend-maxstart+1)/currgroup[0]->len();


				prevgroup[0]=currgroup[0];
				currgroup[0]=currgroup[0]->next_gr;
			}

			float pos_prop=0;
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end + bundledist) { // overlaps positive strand group
				//fprintf(stderr,"\tovlp to pos group: %u-%u\n",prevgroup[2]->start,prevgroup[2]->end);
				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > prevgroup[2]->start ? currgroup[nextgr]->start : prevgroup[2]->start;
					uint minend = currgroup[nextgr]->end < prevgroup[2]->end ? currgroup[nextgr]->end : prevgroup[2]->end;
					if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
					pos_prop+=prevgroup[2]->cov_sum*(minend-maxstart+1)/prevgroup[2]->len();
				}
			}

			while(currgroup[2]!=NULL && currgroup[nextgr]->start <= currgroup[2]->end +bundledist && currgroup[2]->start <= currgroup[nextgr]->end + bundledist) { // overlaps positive strand group
				//fprintf(stderr,"\tovlp to pos group: %u-%u\n",currgroup[2]->start,currgroup[2]->end);

				int grcol = currgroup[2]->color;    // set smallest color for currgroup[$nextgr]
				while(equalcolor[grcol]!=grcol) {
					grcol=equalcolor[grcol];
				}
				currgroup[2]->color=grcol;

				set_strandcol(currgroup[nextgr],currgroup[2],currgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > currgroup[2]->start ? currgroup[nextgr]->start : currgroup[2]->start;
					uint minend = currgroup[nextgr]->end < currgroup[2]->end ? currgroup[nextgr]->end : currgroup[2]->end;
					if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
					pos_prop+=currgroup[2]->cov_sum*(minend-maxstart+1)/currgroup[2]->len();
				}

				prevgroup[2]=currgroup[2];
				currgroup[2]=currgroup[2]->next_gr;
			}

			if(pos_prop) {
				currgroup[nextgr]->neg_prop/=(currgroup[nextgr]->neg_prop+pos_prop);
			}
			else if(currgroup[nextgr]->neg_prop) currgroup[nextgr]->neg_prop=1;
			//fprintf(stderr,"neg_prop=%g pos_prop=%g\n",currgroup[nextgr]->neg_prop,pos_prop);
		}
		else if(nextgr == 0) { // negative strand group
			currgroup[nextgr]->neg_prop=1;
		}
		prevgroup[nextgr]=currgroup[nextgr];
		currgroup[nextgr]=currgroup[nextgr]->next_gr;
    }

	/*
    { // DEBUG ONLY
    	fprintf(stderr,"Colors assigned!\n");
    	for(int sno=0;sno<3;sno++) {
    		fprintf(stderr, "Colors of groups on strand %d:\n",sno);
    		CGroup *procgroup=startgroup[sno];
    		while(procgroup!=NULL) {

    			int grcol = procgroup->color;

    			while(equalcolor[grcol]!=grcol) {
    				grcol=equalcolor[grcol];
    			}

    			int negcol=eqnegcol[grcol];
    			if(eqnegcol[grcol]!=-1){
    				while(equalcolor[negcol]!=negcol) {
    					negcol=equalcolor[negcol];
    				}
    			}
    			int poscol=eqposcol[grcol];
    			if(eqposcol[grcol]!=-1){
    				while(equalcolor[poscol]!=poscol) {
    					poscol=equalcolor[poscol];
    				}
    			}

    			fprintf(stderr, " gr %d(%d,%d,%d,%.6f): %d-%d\n",procgroup->grid,grcol,negcol,poscol,procgroup->cov_sum,procgroup->start,procgroup->end);
    			procgroup=procgroup->next_gr;
    		}
    		//fprintf(stderr,"\n");
    	}
    	//exit(0);
    }
    */



	/*****************************
	 ** Step 8: 'create_bundle' function. 
	 ** 	create bundles : bundles collect connected groups (with same color)
	 *****************************/
	for (int i=0;i<3;i++) {
		currgroup[i]=startgroup[i];
		prevgroup[i]=NULL;
	}

	GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
	GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle : this might be the key for overalps

	GVec<int> group2bundle[3]; // to retrace reads from group no to bundle
	for(int sno=0;sno<3;sno++) {
		group2bundle[sno].Resize(group.Count(),-1);  // for a given group id we get a certain bnode id
		bnode[sno].setFreeItem(false);
	}

	GVec<int> bundlecol(true); // associates a bundle number to a group color
	bundlecol.Resize(equalcolor.Count(),-1);

	while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL) { // there are still groups to process

		int nextgr=get_min_start(currgroup);  // next group based on starting position

		// get group color; I need to redo this to ensure I equalize all colors -> they could still be hanged by set_strandcol
		int grcol = currgroup[nextgr]->color;

		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
		}
		currgroup[nextgr]->color=grcol;

		if(nextgr == 0 || nextgr ==2 || (nextgr==1 &&(eqnegcol[grcol]==-1) && (eqposcol[grcol]==-1))) { // negative or positive strand bundle or unstranded bundle

			int bno=bundlecol[grcol];

			if(bno>-1) { // bundle for group has been created before
				//fprintf(stderr,"Add group=%d to bundle[%d][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
				add_group_to_bundle(currgroup[nextgr],bundle[nextgr][bno],bnode[nextgr],bundledist);
			}
			else { // create new bundle
				bno=create_bundle(bundle[nextgr],currgroup[nextgr],bnode[nextgr]);
				//fprintf(stderr,"Add group=%d to new bundle[%d][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
				bundlecol[grcol]=bno;
			}

			group2bundle[nextgr][currgroup[nextgr]->grid]=bundle[nextgr][bno]->lastnodeid;

		}
		else { // unknown strand : here is where I should compute positive and negative proportions

			if(eqnegcol[grcol]!=-1){
				int negcol=eqnegcol[grcol];
				while(equalcolor[negcol]!=negcol) {
					negcol=equalcolor[negcol];
				}

				int bno=bundlecol[negcol];
				if(bno>-1) { // bundle for group has been created before
					//fprintf(stderr,"Add group=%d to bundle[%d:0][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
					add_group_to_bundle(currgroup[nextgr],bundle[0][bno],bnode[0],bundledist); // this needs bundledist
				}
				else { // create new bundle
					bno=create_bundle(bundle[0],currgroup[nextgr],bnode[0]);
					//fprintf(stderr,"Add group=%d to new bundle[%d:0][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
					bundlecol[negcol]=bno;
				}
				group2bundle[0][currgroup[nextgr]->grid]=bundle[0][bno]->lastnodeid;
			} // if(eqnegcol[grcol]!=-1)

			if(eqposcol[grcol]!=-1){
				int poscol=eqposcol[grcol];
				while(equalcolor[poscol]!=poscol) {
					poscol=equalcolor[poscol];
				}

				int bno=bundlecol[poscol];
				if(bno>-1) { // bundle for group has been created before
					//fprintf(stderr,"Add group=%d to new bundle[%d][%d:2]\n",currgroup[nextgr]->grid,nextgr,bno);
					add_group_to_bundle(currgroup[nextgr],bundle[2][bno],bnode[2],bundledist);
				}
				else { // create new bundle
					bno=create_bundle(bundle[2],currgroup[nextgr],bnode[2]);
					//fprintf(stderr,"Add group=%d to new bundle[%d:2][%d]\n",currgroup[nextgr]->grid,nextgr,bno);
					bundlecol[poscol]=bno;
				}
				group2bundle[2][currgroup[nextgr]->grid]=bundle[2][bno]->lastnodeid;
			}
		}

		currgroup[nextgr]=currgroup[nextgr]->next_gr;

	} // while(currgroup[0]!=NULL || currgroup[1]!=NULL || currgroup[2]!=NULL)


	/*****************************
	 ** Step 9: Clean up no longer needed variables
	 ** 	group.Clear(); maybe I still need this?
	 *****************************/
	equalcolor.Clear();
	eqposcol.Clear();
	eqnegcol.Clear();
	bundlecol.Clear();

/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(after bundles created):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	/*****************************
	 ** Step 10: 'get_covered' function
	 ** 	next variables are in order to remember if I get guide coverage
	 *****************************/
	GVec<int>* bnodeguides=NULL;
	if(bundle[1].Count() && bnode[1].Count()) {
		bnodeguides = new GVec<int>[bnode[1].Count()];
	}

	//if(guides.Count()) fprintf(stderr,"No of guides=%d\n",guides.Count());

	if(c_out || (bundle[1].Count() && bnode[1].Count())) // coverage is needed
		for(int g=0;g<guides.Count();g++) {
			//fprintf(stderr,"consider guide %d\n",g);
			int s=0;
			if(guides[g]->strand=='+') s=2;
			if((c_out && !get_covered(guides[g],bundle[s],bnode[s],junction,NULL,0)) ||
					(bundle[1].Count() && bnode[1].Count() && guides[g]->exons.Count()==1))
				get_covered(guides[g],bundle[1],bnode[1],junction,bnodeguides,g);
		}

	// /*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr, "There are %d unstranded bundles %d negative bundles and %d positive bundles\n",bundle[1].Count(),bundle[0].Count(),bundle[2].Count());
		for(int sno=0;sno<3;sno++) {
			fprintf(stderr, "Bundles on strand %d:\n",sno);
			for(int b=0;b<bundle[sno].Count();b++) {
				int elen=bnode[sno][bundle[sno][b]->lastnodeid]->end-bnode[sno][bundle[sno][b]->startnode]->start+1;
				CBundlenode *currbnode=bnode[sno][bundle[sno][b]->startnode];
				fprintf(stderr,"***Bundle %d with len=%d:",b,elen);
				int nodes=0;
				while(currbnode!=NULL) {
					fprintf(stderr, " %d-%d cov=%.6f",currbnode->start,currbnode->end,currbnode->cov/(currbnode->end-currbnode->start+1));
					currbnode=currbnode->nextnode;
					nodes++;
				}
				currbnode=bnode[sno][bundle[sno][b]->lastnodeid];
				nodes++;
				fprintf(stderr," last node:%d-%d total nodes=%d\n",currbnode->start,currbnode->end,nodes);
			}
		}
	}
	// */

	int geneno=0;

	/*****************************
	 ** Step 11: 'CPrediction': constructor
	 **		predict transcripts for unstranded bundles here
	 *****************************/
	//if(fraglen)
	int g=0;
	for(int b=0;b<bundle[1].Count();b++) { // these are neutral bundles that do not overlap any signed reads

		// I need to address features here too -> TODO

		bool guide_ovlp=false;
		while(g<guides.Count() && (guides[g]->exons.Count()>1 || guides[g]->end<bnode[1][bundle[1][b]->startnode]->start)) {
			g++;
		}
		// now guides[g]->end>=bnode[1][bundle[1][b]->startnode]->start
		if(g<guides.Count() && guides[g]->start<=bnode[1][bundle[1][b]->startnode]->end) guide_ovlp=true;

		if(bundle[1][b]->cov && ((bundle[1][b]->multi/bundle[1][b]->cov)<=mcov*(1-ERROR_PERC) || guide_ovlp || rawreads)) { // && (guides.Count() || adaptive || bundle[1][b]->len >= mintranscriptlen)) { // there might be small transfrags that are worth showing, but here I am ignoring them
    		// bundle might contain multiple fragments of a transcript but since we don't know the complete structure -> print only the pieces that are well represented
    		CBundlenode *currbnode=bnode[1][bundle[1][b]->startnode];
    		int t=1;
    		while(currbnode!=NULL) {
    			//int len=currbnode->end-currbnode->start+1;
    			//float cov=currbnode->cov/(currbnode->end-currbnode->start+1);

    			bool printguides=false;

    			if(!rawreads) for(int i=0;i<bnodeguides[currbnode->bid].Count();i++) {
    				int g=bnodeguides[currbnode->bid][i];
    				geneno++;
    				int glen=guides[g]->end-guides[g]->start+1;
    				if(glen && guides[g]->exons.Count()==1) {
    					RC_TData* tdata=(RC_TData*)(guides[g]->uptr);
    					tdata->in_bundle=3;
    					float gcov=(tdata->t_exons[0])->movlcount/glen;
    					// if(cov<gcov) gcov=cov; WHY DO I DO THIS?? CHECK!!!
    					CPrediction *p=new CPrediction(geneno-1, guides[g], guides[g]->start, guides[g]->end, gcov, guides[g]->strand, glen);
    					if(c_out) {
    						GStr guidecov;
    						guidecov.appendfmt("%.2f",gcov);
    						guides[g]->addAttr("coverage",guidecov.chars());
    						guides[g]->printTranscriptGff(c_out);
    					}
    					GSeg exon(guides[g]->start, guides[g]->end);
    					p->exons.Add(exon);
    					p->exoncov.Add(gcov);
    					if(longreads) p->tlen=-p->tlen;
    					pred.Add(p);
    					printguides=true;
    					guidepred[g]=pred.Count()-1;
    				}
    			}

    			if(!printguides) { // && (adaptive || (cov>=readthr && len>=mintranscriptlen))) {
    				if(t==1) { geneno++;}
    				char sign='.';

    				GVec<CTrimPoint> trimpoint;
    				find_all_trims(refstart,0,currbnode->start, currbnode->end,bpcov,trimpoint); // sign should not matter as I am in a totally neutral zone

    				uint predstart=currbnode->start;
    				uint predend=currbnode->end;

    				for(int i=0;i<trimpoint.Count();i++) {
    					if(trimpoint[i].pos) {
    						if(trimpoint[i].start) {
    							int len=trimpoint[i].pos-CHI_WIN-predstart+1;
    							if(len>mintranscriptlen) {
    								float cov=get_cov(1,predstart-refstart,trimpoint[i].pos-CHI_WIN-refstart,bpcov)/len;

    								//fprintf(stderr,"Store single prediction:%d - %d with cov=%f\n",predstart, trimpoint[i].pos-CHI_WIN, cov);

    								CPrediction *p=new CPrediction(geneno-1, NULL, predstart, trimpoint[i].pos-CHI_WIN, cov, sign, len);
    								GSeg exon(predstart, trimpoint[i].pos-CHI_WIN);

    								p->exons.Add(exon);
    								if(!rawreads) p->exoncov.Add(cov);
    								if(longreads) p->tlen=-p->tlen;
    								pred.Add(p);
    								t++;
    							}
    							predstart=trimpoint[i].pos;
    						}
    						else {
    							int len=trimpoint[i].pos-predstart+1;
    							if(len>mintranscriptlen) {
    								float cov=get_cov(1,predstart-refstart,trimpoint[i].pos-refstart,bpcov)/len;

    			    				//fprintf(stderr,"Store single prediction:%d - %d with cov=%f\n",predstart, trimpoint[i].pos, cov);

    								CPrediction *p=new CPrediction(geneno-1, NULL, predstart, trimpoint[i].pos, cov, sign, len);
    								GSeg exon(predstart, trimpoint[i].pos);

    								p->exons.Add(exon);
    								if(!rawreads) p->exoncov.Add(cov);
    								if(longreads) p->tlen=-p->tlen;
    								pred.Add(p);
    								t++;
    							}
    							predstart=trimpoint[i].pos+CHI_WIN;

    						}
    					}
    				}

    				int len=predend-predstart+1;
    				if(len>mintranscriptlen) {
    					float cov=get_cov(1,predstart-refstart,predend-refstart,bpcov)/len;

    					//fprintf(stderr,"Store single prediction:%d - %d with cov=%f\n",predstart, predend, cov);

    					CPrediction *p=new CPrediction(geneno-1, NULL, predstart, predend, cov, sign, len);
    					GSeg exon(predstart, predend);

    					p->exons.Add(exon);
    					if(!rawreads) p->exoncov.Add(cov);
    					if(longreads) p->tlen=-p->tlen;
    					pred.Add(p);
    					t++;
    				}
    			}
    			currbnode=currbnode->nextnode;
    		}
    	}
    }
    //fprintf(stderr,"Done with unstranded bundles\n");
    if (bnodeguides) delete[] bnodeguides;

	/*****************************
	 ** Step 12: Defining parameters here!!!!
	 ** 	build graphs for stranded bundles here
	 *****************************/
    if(startgroup[0]!=NULL || startgroup[2]!=NULL) {  // Condition 1: there are stranded groups to process

    	// I don't need the groups here anymore : I only use their numbers
    	// group.Clear(); // I will need the proportions later on


		// Create here!!
    	GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
		// Create here!!
    	GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
		// Create here!!
    	CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
		// Create here!!
    	GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
		// Create here!!
    	GVec<int> lastgpos[2];


		// Input variable!!
    	GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
		// Input variable!!
    	GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    	// GVec<int> trnumber[2]; // how many transfrags are on a strand s, in a graph g -> I can find out this from transfrag[s][g].Count()
    	// int ngraph[2]={0,0};   // how many graphs are in each strand: negative (0), or positive(1) -> keep one for each bundle


		// Input variable!!
    	GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
		

    	int bno[2]={0,0};

		/*****************************
		 ** 1. build graph structure
	 	 **     'create_graph_multi', 'construct_treepat_multi'
		 *****************************/
    	for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions

        	// guides appear to be sorted by start --> CHECK THIS!!
        	int g=0;
        	int ng=guides.Count();

    		int s=sno/2; // adjusted strand due to ignoring neutral strand

			int graph_no = bundle[sno].Count();
			// uni_splice_graphGp -> get_gpSize()[s];
			// bundle[sno].Count();

    		char strnd='-';
    		if(s) strnd='+';

    		bundle2graph[s]=NULL;
    		if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    		transfrag[s]=NULL;

    		no2gnode[s] = NULL;

    		tr2no[s]=NULL;
    		gpos[s]=NULL;


			// int refstart = uni_splice_graphGp -> get_refstart();
			// int refend = uni_splice_graphGp -> get_refend();
			// int gpSize = uni_splice_graphGp -> get_gpSize()[s];
			// GVec<int> graphnoGp = uni_splice_graphGp -> get_graphnoGp()[s];
			// GVec<int> edgenoGp = uni_splice_graphGp -> get_graphnoGp()[s];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
			// GPVec<CGraphnode>* no2gnodeGp = uni_splice_graphGp -> get_no2gnodeGp()[s]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i

			// gpSize[s];
			// graphnoGp[s];
			// edgenoGp[s];
			// no2gnodeGp[s];

    		if(graph_no) {
    			transfrag[s]=new GPVec<CTransfrag>[graph_no]; // for each bundle I have a graph ? only if I don't ignore the short bundles


				/****************
				 **  KH Adding comment out
				 ****************/
    			no2gnode[s]=new GPVec<CGraphnode>[graph_no];
				// no2gnode[s] = no2gnodeGp;
				/****************
				 **  END of KH comment out
				 ****************/



    			gpos[s]=new GIntHash<int>[graph_no];


    			GCALLOC(tr2no[s],graph_no*sizeof(CTreePat *));
    			bno[s]=graph_no;

    			for(int b=0;b<graph_no;b++) {



					/****************
					 **  KH Adding comment out
					 ****************/
    				graphno[s].cAdd(0);
					// graphno = uni_splice_graphGp -> get_graphnoGp()[s];
    				edgeno[s].cAdd(0);
					// edgeno = uni_splice_graphGp -> get_edgenoGp()[s];
					/****************
					 **  END of KH comment out
					 ****************/



    				lastgpos[s].cAdd(0);
    				// I am overestmating the edgeno below, hopefully not by too much

    				//fprintf(stderr,"Bundle is: %d - %d start at g=%d sno=%d b=%d\n",bnode[sno][bundle[sno][b]->startnode]->start,bnode[sno][bundle[sno][b]->lastnodeid]->end,g,sno,b);

    				while(g<ng && guides[g]->end<bnode[sno][bundle[sno][b]->startnode]->start) g++;

    				int cg=g;
    				int nolap=0;

					/****************
					 **  This is reference guide. Ignore first.
					 ****************/
    				while(cg<ng && guides[cg]->start<=bnode[sno][bundle[sno][b]->lastnodeid]->end) { // this are potential guides that might overlap the current bundle, and they might introduce extra edges

    					//fprintf(stderr,"...consider guide cg=%d with strand=%c and in_bundle=%d\n",cg,guides[cg]->strand,((RC_TData*)(guides[cg]->uptr))->in_bundle);
    					if((guides[cg]->strand==strnd || guides[cg]->strand=='.') && ((RC_TData*)(guides[cg]->uptr))->in_bundle>=2) {
    						//fprintf(stderr,"Add guide g=%d with start=%d end=%d\n",cg,guides[cg]->start,guides[cg]->end);


							/****************
							 **  KH Adding comment out
							 ****************/
    						edgeno[s][b]+=2; // this is an overestimate: possibly I have both an extra source and an extra sink link
							/****************
							 **  END of KH comment out
							 ****************/


    						nolap++;
    					}
    					cg++;
    				}
					/****************
					 **  This is reference guide. Ignore first.
					 ****************/


    				/*
    				{ // DEBUG ONLY
        				fprintf(stderr,"edgeno[%d][%d]=%d\n",s,b,edgeno[s][b]);
    					if(bundle[sno][b]->cov) {
    						fprintf(stderr,"proc bundle[%d][%d] %f/%f is %f len=%d and %d guides\n",sno,b,
    								bundle[sno][b]->multi,bundle[sno][b]->cov,(float)bundle[sno][b]->multi/bundle[sno][b]->cov,
    								bundle[sno][b]->len,nolap);
    					}
    				}
    				*/

    				// here I can add something in stringtie to lower the mintranscript len if there are guides?
					/*****************************
					 ** bundle is worth processing: it might be that there are small transfrags from source to sink that are worth processing
					 *****************************/
    				if(bundle[sno][b]->cov &&
    						(((bundle[sno][b]->multi/bundle[sno][b]->cov)<=mcov && bundle[sno][b]->len >= mintranscriptlen)
    								||nolap)) { // bundle is worth processing: it might be that there are small transfrags from source to sink that are worth processing

    					/*
    					// first identify drops/spikes in coverage points
    					GVec<CTrimPoint> trims;
    					if(trim) {
    						get_trims(trims,bnode[sno][bundle[sno][b]->startnode],refstart,bpcov);
    						{ // DEBUG ONLY
    							fprintf(stderr,"Trims found:");
    							for(int z=0;z<trims.Count();z++) fprintf(stderr," %d(%d,%f)",trims[z].pos,trims[z].start,trims[z].abundance);
    							fprintf(stderr,"\n");
    						}
    					}
    					graphno[s][b]=create_graph_multi(s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,trims); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure
    					*
    					*/


						/****************
						 **  KH modify
						 ****************/
    					// create graph then
    					graphno[s][b]=create_graph_multi(refstart,s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,gpos,bdata,edgeno[s][b],lastgpos[s][b],guideedge, uni_splice_graphGp); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

// int refstart = uni_splice_graphGp -> get_refstart();
// int refend = uni_splice_graphGp -> get_refend();
// int gpSize = uni_splice_graphGp -> get_gpSize()[s];
// GVec<int> graphnoGp = uni_splice_graphGp -> get_graphnoGp()[s];
// GVec<int> edgenoGp = uni_splice_graphGp -> get_graphnoGp()[s];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
// GPVec<CGraphnode>* no2gnodeGp = uni_splice_graphGp -> get_no2gnodeGp()[s]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i

    					if(graphno[s][b]) tr2no[s][b]=construct_treepat_multi(graphno[s][b],gpos[s][b],transfrag[s][b]);
						/****************
						 **  End of KH modify
						 ****************/



    					else tr2no[s][b]=NULL;
    				}
    				else tr2no[s][b]=NULL;
    			}
    		}
    	}
    	fprintf(stderr,"Done creating graphs\n");
    	// /*
    	{ // DEBUG ONLY
    		printTime(stderr);
    		for(int s=0;s<2;s++) {
    			fprintf(stderr, "There are %d stranded[%d] graphs\n",bno[s],int(2*s));
    			for(int b=0;b<bno[s];b++) {
    				if(graphno[s][b]) {
    					GStr pat;
    					fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges with lastgpos=%d:",int(2*s),b,graphno[s][b],edgeno[s][b],lastgpos[s][b]);
    					for(int nd=1;nd<graphno[s][b]-1;nd++)
    						fprintf(stderr," %d(%d-%d)",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end);
    					fprintf(stderr,"\n");
    					print_pattern(tr2no[s][b],pat,graphno[s][b]);
    				}
    			}
    		}
			fprintf(stderr,"\n");
    	}
    	// */
/*
#ifdef GMEMTRACE
    	double vm,rsm;
    	get_mem_usage(vm, rsm);
    	GMessage("\t\tM(after graphs created):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

		// Move the data cleaning after writing out the DOT file.

		/*****************************
		 ** 2. compute probabilities for stranded bundles
		 **    'get_fragment_pattern' function
		 **        because of this going throu
		 *****************************/
    	for (int n=0;n<readlist.Count();n++) {
	  	/*if(readlist[n]->unitig) { // super-reads are unpaired
    			float srcov=0;
    			for(int i=0;i<readlist[n]->segs.Count();i++)
    				srcov+=get_cov(1,readlist[n]->segs[i].start-refstart,readlist[n]->segs[i].end-refstart,bpcov)/readlist[n]->segs[i].len();
    			if(srcov) get_fragment_pattern(readlist,n,-1,readlist[n]->read_count/srcov,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);

    		}
    		else {*/
    			float single_count=readlist[n]->read_count;
    			for(int j=0; j<readlist[n]->pair_idx.Count();j++) {
    				int np=readlist[n]->pair_idx[j];
    				if(np>-1) {
    					single_count-=readlist[n]->pair_count[j];
    					if(n<np) {
    						get_fragment_pattern(readlist,n,np,readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
    					}
    				}
    			}
    			if(single_count>epsilon) {
    				get_fragment_pattern(readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);
    			}
				// Get the read count for each node 
			//}
    	}


		/*****************************
		 ** 3. Write out global splice graph in DOT format
		*****************************/
		/****************
		 **  KH Adding 
		 ****************/
		if (universal_splice_graph) {
			//  DOT file outut here 
			//  not capacity and rate 
			//  only edge weight
			for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
				int s=sno/2; // adjusted strand due to ignoring neutral strand
				int g_idx = 0;
				for(int b=0;b<bundle[sno].Count();b++) {
					fprintf(stderr, "New writing out place!! Start writing out DOT file!!\n");
					fprintf(stderr,"after traverse:\n");
					if(graphno[s][b]) {
						fprintf(uinigraph_out,"strict digraph %d_%d_%d_%d {", refstart, refend, s, g_idx);
						// graphno[s][b]: number of nodes in graph.
						if(graphno[s][b]) {
							for(int nd=1;nd<graphno[s][b]-1;nd++)
								fprintf(uinigraph_out,"%d[start=%d end=%d cov=%f];",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end,no2gnode[s][b][nd]->cov);

							for(int nd=0;nd<graphno[s][b];nd++) {
								// fprintf(stderr,"Node %d with parents:",i);
								for(int c=0;c<no2gnode[s][b][nd]->child.Count();c++) {
									fprintf(uinigraph_out,"%d->",nd);			
									fprintf(uinigraph_out,"%d;",no2gnode[s][b][nd]->child[c]);
								}
							}
						}
						fprintf(uinigraph_out,"}\n");
						g_idx += 1;
						fprintf(stderr,"g_idx: %d\n", g_idx);
					}
				}
			}
		}
		/****************
		 **  END KH Adding 
		****************/


		/*****************************
		 ** 4. I can clean up some data here:
		*****************************/
		uni_splice_graphGp->Clear();
    	for(int sno=0;sno<3;sno++) {
    		int n=bnode[sno].Count();
    		for(int b=0;b<n;b++) delete bnode[sno][b];
    		bnode[sno].Clear();
    		bundle[sno].Clear();
    	}

/*
#ifdef GMEMTRACE
    	//double vm,rsm;
    	get_mem_usage(vm, rsm);
    	GMessage("\t\tM(read patterns counted):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
    	// shouldn't readlist be also cleared up here? maybe bpcov too?

    	// don't forget to clean up the allocated data here
    	delete [] readgroup;
    	group.Clear();


		/*****************************
		 ** 5. parse graph
		 **    'process_refguides' & 'process_transfrags' & 'find_transcripts' & 'free_treepat'
		 *****************************/
    	for(int s=0;s<2;s++) {

    		for(int b=0;b<bno[s];b++) {
    			// fprintf(stderr,"Process graph[%d][%d] with %d nodes\n",s,b,graphno[s][b]);
    			if(graphno[s][b]) {

    				// include source to guide starts links
    				GVec<CGuide> guidetrf;

    				/*
    				{ // DEBUG ONLY
    					fprintf(stderr,"process refguides for s=%d b=%d edgeno=%d gno=%d lastgpos=%d guidescount=%d\n",s,b,edgeno[s][b],graphno[s][b],lastgpos[s][b],guides.Count());
    					fprintf(stderr,"There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
    					for(int i=0;i<graphno[s][b];i++) {
    						fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len());
    						fprintf(stderr," parents:");
    						for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
    						fprintf(stderr," trf=");
    						for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
    						fprintf(stderr,"\n");
    					}
    				}
    				*/

    				if(guides.Count()) process_refguides(graphno[s][b],edgeno[s][b],gpos[s][b],lastgpos[s][b],no2gnode[s][b],transfrag[s][b],s,guidetrf,bdata);

    				GVec<int> trflong; // non-redundant long transfrags that I can use to guide the long read assemblies
    				//process transfrags to eliminate noise, and set compatibilities, and node memberships
    				process_transfrags(s,graphno[s][b],edgeno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],gpos[s][b],guidetrf,pred,trflong);
    				//get_trf_long(graphno[s][b],edgeno[s][b], gpos[s][b],no2gnode[s][b],transfrag[s][b],geneno,s,pred,trflong);




    				/*
    				{ //DEBUG ONLY
    					//printTime(stderr);
    					fprintf(stderr,"There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
    					for(int i=0;i<graphno[s][b];i++) {
    						fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f capacity=%f rate=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len(), no2gnode[s][b][i]->capacity, no2gnode[s][b][i]->rate);
    						fprintf(stderr," parents:");
    						for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
    						fprintf(stderr," children:");
    						for(int j=0;j<no2gnode[s][b][i]->child.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->child[j]);
    						fprintf(stderr," trf=");
    						for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
    						fprintf(stderr,"\n");
    					}
    					fprintf(stderr,"There are %d transfrags[%d][%d]:\n",transfrag[s][b].Count(),s,b);
    					for(int t=0;t<transfrag[s][b].Count();t++) {
    						fprintf(stderr,"%d: ",t);
    						//printBitVec(transfrag[s][b][t]->pattern);
    						fprintf(stderr," %f(%f) long=%d short=%d nodes=%d",transfrag[s][b][t]->abundance,transfrag[s][b][t]->srabund, transfrag[s][b][t]->longread,transfrag[s][b][t]->shortread,transfrag[s][b][t]->nodes.Count());
    						for(int i=0;i<transfrag[s][b][t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[s][b][t]->nodes[i]);
    						if(!transfrag[s][b][t]->abundance) fprintf(stderr," *");
    						fprintf(stderr,"\n");
    					}

    				}
    				*/

/*
#ifdef GMEMTRACE
    				double vm,rsm;
    				get_mem_usage(vm, rsm);
    				GMessage("\t\tM(after process_transfrags):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

    				// fprintf(stderr,"guidetrf no=%d\n",guidetrf.Count());

    				//if(!longreads) {
    				// find transcripts now
    				if(!rawreads) geneno=find_transcripts(graphno[s][b],edgeno[s][b],gpos[s][b],no2gnode[s][b],transfrag[s][b],
    						geneno,s,guidetrf,guides,guidepred,bdata,trflong);
    				//}
    				for(int g=0;g<guidetrf.Count();g++) delete guidetrf[g].trf;

    				/*
    				{ //DEBUG ONLY
    					printTime(stderr);
    					fprintf(stderr,"Processed transcripts for s=%d b=%d\n",s,b);
    				}
    				*/

/*
#ifdef GMEMTRACE
    				//double vm,rsm;
    				get_mem_usage(vm, rsm);
    				GMessage("\t\tM(after processed transcripts):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
    			}
    			// clean up what can be cleaned
    			if(tr2no[s][b]) free_treepat(tr2no[s][b]);
    		}

    		// final clean up: no2gnode, no2tr, transfrag, bundle2graph
    		if(bundle2graph[s]) delete [] bundle2graph[s];
    		if(transfrag[s]) delete [] transfrag[s];
    		if(no2gnode[s]) delete [] no2gnode[s];
    		if(gpos[s]) delete [] gpos[s];
    		if(tr2no[s]) GFREE(tr2no[s]);
    	}

	} // end if(startgroup[0]!=NULL || startgroup[2]!=NULL)
    else { 	// Condition 2: no stranded groups to process

    	delete [] readgroup;
    	// clean up readgroup, bundle
    	for(int sno=0;sno<3;sno++) {
    		int n=bnode[sno].Count();
    		for(int b=0;b<n;b++) delete bnode[sno][b];
    		bnode[sno].Clear();
    	}
    }

/*
#ifdef GMEMTRACE
    //double vm,rsm;
    get_mem_usage(vm, rsm);
	GMessage("\t\tM(e):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

    /*
    { // DEBUG ONLY
    	for(int i=0;i<pred.Count();i++) {
    		if(pred[i]->t_eq) fprintf(stderr,"%s ",pred[i]->t_eq->getID());
    		fprintf(stderr,"pred[%d] (cov=%f,strand=%c):",i,pred[i]->cov,pred[i]->strand);
    		for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
    		fprintf(stderr,"\n");
    	}
    }
    */

    // don't forget to clean up the allocated data here
    return(geneno);
}


int build_graphs_unispg(BundleData* bdata, UniSpliceGraphGp* uni_splice_graphGp) {
	int refstart = bdata->start;
	int refend = bdata->end+1;
	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GPVec<GffObj>& guides = bdata->keepguides;
	GVec<float>* bpcov = bdata->bpcov; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles
	GList<CPrediction>& pred = bdata->pred;



	/****************
	 **  These are parameters initialized to build graphs. 
	 ****************/
	// form groups on strands: all groups below are like this: 0 = negative strand; 1 = unknown strand; 2 = positive strand
	GPVec<CGroup> group;
	CGroup *currgroup[3]={NULL,NULL,NULL}; // current group of each type
	CGroup *startgroup[3]={NULL,NULL,NULL}; // start group of each type
	int color=0; // next color to assign
	GVec<int> merge; // remembers merged groups
	GVec<int> equalcolor; // remembers colors for the same bundle
	GVec<int> *readgroup=new GVec<int>[readlist.Count()]; // remebers groups for each read; don't forget to delete it when no longer needed
	// parameter -e ; for mergeMode includes estimated coverage sum in the merged transcripts
	GVec<int> guidepred; // for eonly keeps the prediction number associated with a guide
	GArray<GEdge> guideedge; // 0: negative starts; 1 positive starts
	/*
	GPVec<GPtFeature>& feature = bdata->ptfs; // these are point features (confirmed starts/stops)

	for(int i=0;i<feature.Count();i++) {
		if(feature[i]->ftype==GPFT_TSS)
			fprintf(stderr,"TSS at position %d on strand %d\n",feature[i]->coord,feature[i]->strand);
		if(feature[i]->ftype==GPFT_CPAS)
			fprintf(stderr,"CPAS at position %d on strand %d\n",feature[i]->coord,feature[i]->strand);
	}
	*/
	/****************
	 **  These are parameters initialized to build graphs. 
	 ****************/
	fprintf(stderr,"build_graphs with %d guides\n",guides.Count());

	/*****************************
	 ** Step 1: this part is for adjusting leftsupport and rightsupport when considering all junctions that start at a given point
	 ** 	sort junctions -> junctions are sorted already according with their start, but not their end
	 *****************************/
	GList<CJunction> ejunction(junction);
	ejunction.setFreeItem(false);
	if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);

	uint start=0;
	uint end=0;
	double leftsupport[2]={0,0};    // should be strand based
	double rightsupport[2]={0,0};
	bool higherr=false;
	char leftcons=-1;
	char rightcons=-1;

	/*****************************
	 ** Step 2: there are some reads that contain very bad junctions -> need to find better closest junctions
	 *****************************/
	if(higherr) { 
		uint juncsupport=junctionsupport;
		if(longreads)
			juncsupport=sserror;
		else if(mixedMode) juncsupport=sserror/DROP;
		//fprintf(stderr,"In higherr!\n");
		GVec<int> jstarts; // good starting junctions
		GVec<int> jends; // good ending junctions
		if(viral) {
			for(int i=1;i<junction.Count();i++) { // junction is sorted based on start

				if(junction[i]->strand && junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					if(junction[i]->nreads_good>=0 && (junction[i]->nreads_good<1.25*junctionthr || !good_junc(*junction[i],refstart,bpcov))) { // threshold for bad junctions is higher; (should I also add that too short junctions not to be accepted?)
						//junction[i]->strand=0; // just delete junction if it's low count
						junction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}

					int j=i-1;
					while(j>0 && junction[i]->start-junction[j]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand && junction[j]->nreads_good>=0 && junction[i]->nreads<junction[j]->nreads &&
								abs((int)junction[i]->end-(int)junction[j]->end)<(int)juncsupport) { // j was not elminated and is better
							if(junction[i]->nreads_good<0) { // i was eliminated before
								int k=-junction[i]->nreads_good;
								if(junction[k]->nreads<junction[j]->nreads) junction[i]->nreads_good=-j;
							}
							else {
								junction[i]->nreads_good=-j;
							}
						}
						j--;
					}
				}
			}
			for(int i=junction.Count()-1;i>0;i--) {
				if(junction[i]->strand && junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					int j=i+1;
					while(j<junction.Count() && junction[j]->start-junction[i]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand && junction[j]->nreads_good>=0 && junction[i]->nreads<junction[j]->nreads &&
								abs((int)junction[i]->end-(int)junction[j]->end)<(int)juncsupport) { // j was not elminated and is better
							if(junction[i]->nreads_good<0) { // i was eliminated before
								int k=-junction[i]->nreads_good;
								if(junction[k]->nreads<junction[j]->nreads) junction[i]->nreads_good=-j;
							}
							else {
								junction[i]->nreads_good=-j;
							}
						}
						j++;
					}
				}
			}
		} else {

			float tolerance=1-ERROR_PERC;

			// strand based version
			for(int i=1;i<junction.Count();i++) { // junction is sorted based on start

				//fprintf(stderr,"junct[%d]:%d-%d:%d lefttsupport=%f nm=%f mm=%f nreads=%f nreads_good=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->leftsupport,junction[i]->nm,junction[i]->mm,junction[i]->nreads,junction[i]->nreads_good);
				if(junction[i]->strand) {
					if(junction[i]->nm && !junction[i]->guide_match && junction[i]->nm>=junction[i]->nreads) { // this is a bad junction -> check if it's maximal;
					if(junction[i]->nreads_good>=0 && junction[i]->nreads_good<1.25*junctionthr) { // threshold for bad junctions is higher; (should I also add that too short junctions not to be accepted?)
						//junction[i]->strand=0; // just delete junction if it's low count
						junction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}

					int j=i-1;
					float support=0;
					bool searchjunc=true;
					bool reliable=false;
					//if(j>=0) fprintf(stderr,"...start at junct:%d-%d:%d leftsupport=%f dist=%d\n",junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport,junction[i]->start-junction[j]->start);
					while(j>0 && junction[i]->start-junction[j]->start<juncsupport) {
						if(junction[j]->strand==junction[i]->strand) {
							if(junction[j]->start==junction[i]->start) { // found a junction with the same start -> I have already searched it if it's bad
								if(junction[j]->nreads<0) {
									junction[i]->nreads=junction[j]->nreads;
									searchjunc=false;
								}
								break;
							}
							else if(junction[j]->guide_match || junction[j]->nm<junction[j]->nreads) { // nearby junction is much more reliable
								if(junction[j]->leftsupport>junction[i]->leftsupport*tolerance) { // the good junction is close enough
									reliable=true;
									junction[i]->nreads=-j;
									support=junction[j]->leftsupport;
									break;
								}
							}
							else if(junction[j]->leftsupport>support && junction[i]->start-junction[j]->start<sserror && junction[j]->leftsupport*tolerance>junction[i]->leftsupport) {
								//fprintf(stderr,"...1 compare to [%d]:%d-%d:%d leftsupport=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport);
								junction[i]->nreads=-j;
								support=junction[j]->leftsupport;
							}
						}
						j--;
					}
					if(searchjunc) {
						j=i+1;
						int dist=juncsupport;
						if(junction[i]->nreads<0) {
							dist=junction[i]->start-junction[abs((int)junction[i]->nreads)]->start;
						}
						//if(j<junction.Count()) fprintf(stderr,"...start at junct:%d-%d:%d leftsupport=%f dist=%d\n",junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport,junction[j]->start-junction[i]->start);
						while(j<junction.Count() && junction[j]->start-junction[i]->start<juncsupport) {
							if(junction[j]->strand==junction[i]->strand && junction[j]->start!=junction[i]->start) {
								int d=(int)(junction[j]->start-junction[i]->start);
								if(junction[j]->guide_match || junction[j]->nm<junction[j]->nreads) {
									if((d<dist || (d==dist && junction[j]->leftsupport>support)) && junction[j]->leftsupport>junction[i]->leftsupport*tolerance) {
										junction[i]->nreads=-j;
										support=junction[j]->leftsupport;
										break;
									}
								}
								else if(!reliable && junction[j]->leftsupport>support && (uint)d<sserror && junction[j]->leftsupport*tolerance>junction[i]->leftsupport) { // junction is not best within window
									//fprintf(stderr,"...2 compare to [%d]:%d-%d:%d leftsupport=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j]->leftsupport);
									junction[i]->nreads=-j;
									support=junction[j]->leftsupport;
								}
							}
							j++;
						}
					}
				}
					else if(mixedMode && junction[i]->nm<junction[i]->nreads) {
						jstarts.Add(i);
					}
				}
				//fprintf(stderr,"ejunct[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);

				if(ejunction[i]->strand) {
					if(ejunction[i]->nm && !ejunction[i]->guide_match && ejunction[i]->nm>=ejunction[i]->nreads) { // this is a bad junction -> check if it's maximal
					if(ejunction[i]->nreads_good>=0 && ejunction[i]->nreads_good<1.25*junctionthr) { // threshold for bad junctions is higher
						//ejunction[i]->strand=0;
						ejunction[i]->mm=-1;
						//fprintf(stderr,"...delete due to being under threshold\n");
					}
					int j=i-1;
					float support=0;
					bool searchjunc=true;
					bool reliable=false;
					//if(j>=0) fprintf(stderr,"...start at junct:%d-%d:%d rightsupport=%f dist=%d\n",ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport,ejunction[i]->end-ejunction[j]->end);
					while(j>0 && ejunction[i]->end-ejunction[j]->end<juncsupport) {
						if(ejunction[j]->strand==ejunction[i]->strand) {
							if(ejunction[j]->end==ejunction[i]->end) {
								if(ejunction[j]->nreads_good<0) {
									ejunction[i]->nreads_good=ejunction[j]->nreads_good;
									searchjunc=false;
								}
								break;
							}
							else if(ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads) { // nearby junction is much more reliable
								if(ejunction[j]->rightsupport>ejunction[i]->rightsupport*tolerance) { // the good junction is close enough
									reliable=true;
									ejunction[i]->nreads_good=-j;
									support=ejunction[j]->rightsupport;
									break;
								}
							}
							else if(ejunction[j]->rightsupport>support && ejunction[i]->end-ejunction[j]->end < sserror && ejunction[j]->rightsupport*tolerance>ejunction[i]->rightsupport) {
								//fprintf(stderr,"...1 compare to [%d]:%d-%d:%d rightsupport=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport);
								ejunction[i]->nreads_good=-j;
								support=ejunction[j]->rightsupport;
							}
						}
						j--;
					}
					if(searchjunc) {
						j=i+1;
						int dist=juncsupport;
						if(ejunction[i]->nreads_good<0) {
							dist=ejunction[i]->end-ejunction[abs((int)ejunction[i]->nreads_good)]->end;
						}
						//if(j<junction.Count()) fprintf(stderr,"...start at junct:%d-%d:%d rightsupport=%f dist=%d\n",ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport,ejunction[j]->end-ejunction[i]->end);
						while(j<junction.Count() && ejunction[j]->end-ejunction[i]->end<juncsupport) {
							if(ejunction[j]->strand==ejunction[i]->strand && ejunction[j]->end!=ejunction[i]->end) {
								int d=ejunction[j]->end-ejunction[i]->end;
								if(ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads) {
									if((d<dist || (d==dist && ejunction[j]->rightsupport>support)) && ejunction[j]->rightsupport>ejunction[i]->rightsupport*tolerance) {
										ejunction[i]->nreads_good=-j;
										support=ejunction[j]->rightsupport;
										break;
									}
								}
								else if((!reliable && ejunction[j]->rightsupport>support && ejunction[j]->end-ejunction[i]->end < sserror && ejunction[j]->rightsupport*tolerance>ejunction[i]->rightsupport) || ((int)(ejunction[j]->end-ejunction[i]->end)<dist && (ejunction[j]->guide_match || ejunction[j]->nm<ejunction[j]->nreads))) {
									//fprintf(stderr,"...2 compare to [%d]:%d-%d:%d rightsupport=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j]->rightsupport);
									ejunction[i]->nreads_good=-j;
									support=ejunction[j]->rightsupport;
								}
							}
							j++;
						}
					}
				}
					else if(mixedMode && ejunction[i]->nm<ejunction[i]->nreads) {
						jends.Add(i);
					}
				}
			}
		}
		if(mixedMode) { // check if there are junctions inside bigger junctions that can form small exons
			int si=0;
			for(int ei=0;ei<jends.Count();ei++) {
				while(si<jstarts.Count() && junction[jstarts[si]]->start<=ejunction[jends[ei]]->end) si++;
				int k=si;
				while(k<jstarts.Count() && junction[jstarts[k]]->start-ejunction[jends[ei]]->end<SMALL_EXON) {
					if(junction[jstarts[k]]->strand == ejunction[jends[ei]]->strand) {
						CJunction jn(ejunction[jends[ei]]->start,junction[jstarts[k]]->end,junction[jstarts[k]]->strand);
						int oidx=-1;
						if (junction.Found(&jn, oidx) && junction[oidx]->nm>=junction[oidx]->nreads) { // candidate junction for deletion
							junction[oidx]->strand=0;
							break;
						}
					}
					k++;
				}
			}
		}
	} //if higherr
	/*
	{ // DEBUG ONLY
		for(int i=0;i<junction.Count();i++) {
			if(junction[i]->guide_match) fprintf(stderr,"G");
			if(junction[i]->strand && junction[i]->nreads>0 && junction[i]->nreads_good>0 && junction[i]->mm>=0) fprintf(stderr,"***");
			//if(junction[i]->strand && junction[i]->mm>=0)
				fprintf(stderr,"Junction[%d] %d-%d:%d has nm=%f mm=%f nreads=%f nreads_good=%f leftsupport=%f and rightsupport=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,
					junction[i]->nm,junction[i]->mm,junction[i]->nreads,
					junction[i]->nreads_good,junction[i]->leftsupport,junction[i]->rightsupport);
		}
	}
	*/
	//int **readgroup = new int*[readlist.Count()];
/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	/*****************************
	 ** Step 3: junctions filtering & sort
	 *****************************/
	//float fraglen=0;
	//uint fragno=0;
	//GHash<bool> boundaryleft;
	//GHash<bool> boundaryright;
	GIntHash<bool> boundaryleft;
	GIntHash<bool> boundaryright;

	bool resort=false;
	int njunc=junction.Count();
	for (int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		// version that does not adjust read but discards it completely
		bool keep=true;
		int i=0;

		// /*
		fprintf(stderr,"check read[%d]:%d-%d:%d refstart=%d w/exons:",n,rd.start,rd.end,rd.strand,refstart);
		for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
		fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);
		i=0;
		// */

		while(i<rd.juncs.Count()) {
			CJunction& jd=*(rd.juncs[i]);
			fprintf(stderr, " read junc %d-%d:%d nreads=%f nreads_good=%f nm=%f support=%f,%f between exons:%d-%d and %d-%d\n", jd.start, jd.end, jd.strand,jd.nreads,jd.nreads_good,jd.nm,jd.leftsupport,jd.rightsupport,rd.segs[i].start,rd.segs[i].end,rd.segs[i+1].start,rd.segs[i+1].end);

			bool changeright=jd.nreads_good<0;
			bool changeleft=jd.nreads<0;
			if(viral && changeright) {
				changeleft=true;
				jd.nreads=jd.nreads_good;
			}

			if(jd.strand) {
				// Definition of `double mm`; // number of reads that support a junction with both anchors bigger than longintronanchor
				if(!good_junc(jd,refstart,bpcov)) jd.mm=-1; // bad junction
				if(changeleft || changeright) jd.strand=0;
				if(jd.mm<0) jd.strand=0; // I can't believe the junction because it's too low
			}

			if(!jd.strand) { // found a bad junction
				//instead of deleting read -> just delete junction and make sure that add_read_to_group adds read to multiple groups
				//jd.nreads-=rd.read_count;

				//if(jd.mm>=0 && (changeleft || changeright)) { // not sure why jd.mm should be positive here since I am searching the junction against the others anyway
				if(changeleft || changeright) {
					uint newstart=rd.segs[i].end;
					uint newend=rd.segs[i+1].start;
					int jk=-1;
					int ek=-1;
					if(changeleft) {
						jk=abs(int(jd.nreads));
						if(junction[jk]->nreads<0) { // not a valid left support -> delete junction
							newstart=rd.segs[i].start-1;
						}
						else {
							newstart=junction[jk]->start; // junction jk is good
							//fprintf(stderr,"junction has newstart=%d from junction[jk=%d]\n",newstart,jk);
							if(junction[jk]->strand) jk=-1;
						}
					}
					if(changeright) {
						ek=abs(int(jd.nreads_good));
						if(viral) {
							if(junction[ek]->nreads_good<0) {
								newend=rd.segs[i+1].end+1;
							}
							else {
								newend=junction[ek]->end;
								//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
							}
						}
						else {
							if(ejunction[ek]->nreads_good<0) {
								newend=rd.segs[i+1].end+1;
							}
							else {
								newend=ejunction[ek]->end;
								//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
								if(ejunction[ek]->strand) ek=-1;
							}
						}
					}
					//if(jd.mm>=0 && newstart>=rd.segs[i].start && newend<=rd.segs[i+1].end) { // junction inside read boundaries
					if(newstart>=rd.segs[i].start && newend<=rd.segs[i+1].end && newstart<=newend) { // junction inside read boundaries
						bool searchjunc=true;
						bool addjunction=true;
						rd.segs[i].end=newstart;   // adjust start
						rd.segs[i+1].start=newend; // adjust end
						if(jd.mm>=0) {
							if(jk>0) { // search junctions
								int k=jk;
								while(k>0 && junction[k]->start==newstart) {
								   if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
								   if(rd.strand==junction[k]->strand && junction[k]->end==newend) {
									rd.juncs.Put(i,junction[k]);
										searchjunc=false;
										break;
									}
									k--;
								}
								if(searchjunc) {
									k=jk+1;
									while(k<junction.Count() && junction[k]->start==newstart) {
										if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
										if(rd.strand==junction[k]->strand && junction[k]->end==newend) {
											rd.juncs.Put(i,junction[k]);
											searchjunc=false;
											break;
										}
										k++;
									}
								}
								if(searchjunc) { // I did not find junction -> I need to create a new one
									// first check if I already added such a junction
									for(k=njunc;k<junction.Count();k++) {
										if(rd.strand==junction[k]->strand && junction[k]->start==newstart && junction[k]->end==newend) {
											rd.juncs.Put(i,junction[k]);
											searchjunc=false;
											break;
										}
									}
									if(searchjunc && addjunction) {
										if(!resort) {
											junction.setSorted(false);
											ejunction.setSorted(false);
											resort=true;
										}
										//fprintf(stderr,"Add new junction:%d-%d:%d at position %d njunc=%d\n",newstart,newend,rd.strand,junction.Count(),njunc);
										CJunction *junc=new CJunction(newstart,newend,rd.strand);
										junction.Add(junc);
										ejunction.Add(junc);
										rd.juncs.Put(i,junc);
									}
								}
							}
							else if(ek>0) { // ek>0 => search ejunctions; ek>0 because I have either changeleft or changeright
								int k=ek;
								while(k>0 && ejunction[k]->end==newend) {
									if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
									if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart) {
										rd.juncs.Put(i,ejunction[k]);
										searchjunc=false;
										break;
									}
									k--;
								}
								if(searchjunc) {
									k=ek+1;
									while(k<ejunction.Count() && ejunction[k]->end==newend) {
										if(mixedMode && junction[k]->nm<junction[k]->nreads) addjunction=false;
										if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart) {
											rd.juncs.Put(i,ejunction[k]);
											searchjunc=false;
											break;
										}
										k++;
									}
								}
								if(searchjunc && addjunction) { // I did not find junction -> I need to create a new one
									// first check if I already added such a junction
									for(k=njunc;k<ejunction.Count();k++) {
										if(rd.strand==ejunction[k]->strand && ejunction[k]->start==newstart && ejunction[k]->end==newend) {
											rd.juncs.Put(i,ejunction[k]);
											searchjunc=false;
											break;
										}
									}
									if(searchjunc) {
										if(!resort) {
											junction.setSorted(false);
											ejunction.setSorted(false);
											resort=true;
										}
										//fprintf(stderr,"Add new junction:%d-%d:%d at position %d njunc=%d\n",newstart,newend,rd.strand,junction.Count(),njunc);
										CJunction *junc=new CJunction(newstart,newend,rd.strand);
										junction.Add(junc);
										ejunction.Add(junc);
										rd.juncs.Put(i,junc);
									}
								}
							}
						}
					}
					//fprintf(stderr, "read[%d] adjusted to junction:%d-%d\n",n,rd.segs[i].end,rd.segs[i+1].start);
				}

				// because read might be poorly mapped I also have to unpair it
				if(keep) { // this is the first time I unpair read -> remove strand if pair has single exon?
					for(int p=0;p<rd.pair_idx.Count();p++) {
						int np=rd.pair_idx[p];
						if(np>-1) {
							rd.pair_idx[p]=-1;
							for(int j=0;j<readlist[np]->pair_idx.Count();j++) // also unpair read np to n
								if(readlist[np]->pair_idx[j]==n) {
									readlist[np]->pair_idx[j]=-1;
									break;
								}
							// remove strand for pair read if single exon ? is strand assigned before here ?
							if(!readlist[np]->juncs.Count()) readlist[np]->strand=0;
						}
					}
					keep=false;
				}
			}
			else if(guides.Count()){ // need to remember boundary
				bool exist=true;
		    	//GStr bs((int)jd.start);
		    	if(!boundaryleft[jd.start]) boundaryleft.Add(jd.start,exist);
		    	//GStr be((int)jd.end);
		    	if(!boundaryright[jd.end]) boundaryright.Add(jd.end,exist);
			}
			i++;

		}

		if(!keep) { // read has bad junctions -> check if I need to unstrand it
			if(rd.strand) {
				bool keepstrand=false;
				for(int j=0;j<rd.juncs.Count();j++) if(rd.juncs[j]->strand) { keepstrand=true; break; }
				if(!keepstrand) rd.strand=0;
			}
		}
		else if(!rd.juncs.Count() && rd.strand) { // check if I need to unstrand current read due to future mapping
			if(rd.pair_idx.Count()) {
				bool keepstrand=false;
				for(int p=0;p<rd.pair_idx.Count();p++) {
					int np=rd.pair_idx[p];
					if(np>-1 && n<np) {
						if(!readlist[np]->juncs.Count()) {
							if(readlist[np]->strand==rd.strand) {
								keepstrand=true;
								break;
							}
						}
						else {
							for(int j=0;j<readlist[np]->juncs.Count();j++) {
								CJunction& jd=*(readlist[np]->juncs[j]);
								if(jd.strand && good_junc(jd,refstart,bpcov)) {
									keepstrand=true;
									break;
								}
							}
							if(!keepstrand) readlist[np]->strand=0;
						}
					}
					if(keepstrand) break;
				}
				if(!keepstrand) rd.strand=0;
			}
		}


		//if(rd.juncs.Count()) fprintf(stderr,"read[%d] keep=%d\n",n,keep);
		//if(rd.strand) fprintf(stderr,"read[%d] has strand %d\n",n,rd.strand);


		//if(keep) { // if it's a good read that needs to be kept


			/*fprintf(stderr,"add read %d:%d-%d w/count=%g for color=%d with npairs=%d\n",n,readlist[n]->start,readlist[n]->end,readlist[n]->read_count,color,readlist[n]->pair_idx.Count());
			fprintf(stderr,"add read[%d]:%d-%d:%d w/count=%g w/exons:",n,readlist[n]->start,readlist[n]->end,readlist[n]->strand,readlist[n]->read_count);
			for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
			fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);*/

			color=add_read_to_group(n,readlist,color,group,currgroup,startgroup,readgroup,equalcolor,merge);

			// count fragments
			if(!rd.unitig)
				bdata->frag_len+=rd.len*rd.read_count; // TODO: adjust this to work with FPKM for super-reads and Pacbio
			double single_count=rd.read_count;
			if(keep) for(int i=0;i<rd.pair_idx.Count();i++) {
				// I am not counting the fragment if I saw the pair before and it wasn't deleted
				if(rd.pair_idx[i]!=-1 && n>rd.pair_idx[i] && readlist[rd.pair_idx[i]]->nh) {// only if read is paired and it comes first in the pair I count the fragments
					single_count-=rd.pair_count[i];
				}
			}
			if(!rd.unitig && single_count>epsilon) {
				bdata->num_fragments+=single_count; // TODO: FPKM will not work for super-reads here because I have multiple fragments in
												    // a super-read -> I might want to re-estimate this from coverage and have some input for read length; or I might only use TPM
			}

			//fprintf(stderr,"now color=%d\n",color);
		//}
		//else { fprintf(stderr,"read[%d] is not kept\n",n);}
		//else clean_read_junctions(readlist[n]);
	}

	if(resort) {
		junction.setSorted(true);
		ejunction.setSorted(juncCmpEnd);
	}
	//fprintf(stderr,"fragno=%d fraglen=%g\n",fragno,fraglen);
	//if(fragno) fraglen/=fragno;

	return 0;
}


int infer_transcripts_multi(BundleData* bundle, UniSpliceGraphGp* uni_splice_graphGp) {
	int geneno=0;

	//DEBUG ONLY: 	showReads(refname, readlist);

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(s):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/
	if(mergeMode) {
	// 	//count_merge_junctions(bundle->readlist,bundle->covflags); // this make sense only if I want to count junctions
		geneno = build_merge(bundle);
	}
	else if(bundle->keepguides.Count() || !eonly) {
		//fprintf(stderr,"Process %d reads from %lu.\n",bundle->readlist.Count(),bundle->numreads);

		count_good_junctions(bundle);

		// geneno = build_graphs_multi(bundle, uni_splice_graphGp);
		geneno = build_graphs_unispg(bundle, uni_splice_graphGp);
	}


/*
#ifdef GMEMTRACE
	//double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(e):infer_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	return(geneno);
}