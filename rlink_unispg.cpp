#include "rlink.h"
#include "rlink_unispg.h"
#include "GBitVec.h"
#include <float.h>
#include <cmath>

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

extern bool multiMode;
extern bool unispgMode;

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

extern GStr label;

inline int edge(int min, int max, int gno) {
	//return((gno-1)*min-min*(min-1)/2+max-min); // this should be changed if source to node edges are also stored
	return((gno-1)*(min+1)-min*(min-1)/2+max-min); // this includes source to node edges
}


CGraphnode *create_graphnode_unispg(int s, int g, uint start,uint end,int nodeno,CBundlenode *bundlenode,
		GVec<CGraphinfo> **bundle2graph,GPVec<CGraphnode> **no2gnode) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"create_graphnode_unispg[%d][%d]:%d-%d nodeno=%d\n",s,g,start,end,nodeno);
	}
	*/
	CGraphnode* gnode=new CGraphnode(start,end,nodeno);
	CGraphinfo ginfo(g,nodeno);
	bundle2graph[s][bundlenode->bid].Add(ginfo);
	no2gnode[s][g].Add(gnode);

	return(gnode);
}

GBitVec traverse_dfs_unispg(int s,int g,CGraphnode *node,CGraphnode *sink,GBitVec parents,int gno, GVec<bool>& visit,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag, int &edgeno,GIntHash<int> **gpos,int &lastgpos){

	//fprintf(stderr,"Traverse node %d\n",node->nodeid);

	if(visit[node->nodeid]) {
		node->parentpat = node->parentpat | parents;
		for(int n=0;n<gno;n++) {
			if(parents[n]) // add node's children to all parents of node
				no2gnode[s][g][n]->childpat = no2gnode[s][g][n]->childpat | node->childpat;
			else if(node->childpat[n])
				no2gnode[s][g][n]->parentpat = no2gnode[s][g][n]->parentpat | node->parentpat;
		}
	}
	else {
		node->childpat.resize(gno+edgeno);
		node->parentpat.resize(gno+edgeno);
		node->parentpat = node->parentpat | parents;
		visit[node->nodeid]=true;
		parents[node->nodeid]=1; // add the node to the parents

		if(node->parent.Count()==1 && !node->parent[0]) { // node has source only as parent -> add transfrag from source to node
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

	    for(int i=0; i< n; i++) { // for all children
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
	    	node->childpat = node->childpat | traverse_dfs_unispg(s,g,no2gnode[s][g][node->child[i]],sink,childparents,gno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	    }
	} // end else from if(visit[node->nodeid])

	GBitVec children = node->childpat;
	children[node->nodeid]=1;

	return(children);
}


int create_graph_unispg(int refstart,int s,int g,CBundle *bundle,GPVec<CBundlenode>& bnode,
		GList<CJunction>& junction,GList<CJunction>& ejunction,GVec<CGraphinfo> **bundle2graph,
		GPVec<CGraphnode> **no2gnode,GPVec<CTransfrag> **transfrag,GIntHash<int> **gpos,BundleData* bdata,
		int &edgeno,int &lastgpos,GArray<GEdge>& guideedge, int refend=0){

	if (multiMode) {
        // int uni_refstart = unispg_gp -> get_refstart();
		// int uni_refend = unispg_gp -> get_refend();
        // fprintf(stderr, "* uni_refstart: %d\n", uni_refstart);
        // fprintf(stderr, "* uni_refend: %d\n", uni_refend);
	} else if (unispgMode){
		// int* uni_gpSize = unispg_gp -> get_gpSize();
		// GVec<int>* uni_graphnoGp = unispg_gp -> get_graphnoGp();
		// GVec<int>* uni_edgenoGp = unispg_gp -> get_edgenoGp();  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
		// GPVec<CGraphnode>** uni_no2gnodeGp = unispg_gp -> get_no2gnodeGp(); // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
	} else {
		// normalMode
	}

	/****************
	 **  KH Adding 
	****************/
	fprintf(stdout, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	fprintf(stdout, "&&&&&& Start 'create_graph_unispg'\n");
	fprintf(stdout, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	/****************
	 **  END KH Adding 
	****************/

	GVec<float>* bpcov = bdata ? bdata->bpcov : NULL; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();
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

	int graphno=1; // number of nodes in graph
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
	 ** Step 3: 'create_graphnode_unispg' function
	 ** 	Process nodes in the bundle.
	 *****************************/
	int f=0; // feature index
	uint bundle_start=bundlenode->start;
	uint bundle_end=bnode[bundle->lastnodeid]->end;
	GHashMap<int, int> global2local_nodehash(false); //hash of pointers
	global2local_nodehash.Add(0, 0);
	// I want to process bundles & compare the local to the global graph.
	int nd_global=1;

	/***************************
	 ** In C++ visualization
	***************************/
	// vector<vector<int> > exonIntervals;
	// vector<vector<int> > intronIntervals;
	// vector<vector<int> > exonIntervals_unispg;
	// vector<int> exon_tmp;
	// vector<int> intron_tmp;
	/***************************
	 ** In C++ visualization
	***************************/

	while(bundlenode!=NULL) {
		fprintf(stderr,"process bundlenode %d-%d:%d bpcov_count=%d refstart=%d\n",bundlenode->start,bundlenode->end,s,bpcov->Count(),refstart);

		uint currentstart=bundlenode->start; // current start is bundlenode's start
		uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now:q
		int end=0;
		while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
			if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
				end=1;
			}
			nje++;
		}

		GVec<CPred> lstart; // CPred: prediction point class
		GVec<CPred> lend;
		int fs=-1; // first start feature index in lstart
		int fe=-1; // first end feature index in lend

		// see if I need to adjust the start to ignore little hanging pieces that make no sense
		if(!end) {
			while(nje<njunctions && ejunction[nje]->strand+1!=2*s) nje++; // skip all junctions that are not on the same strand
			if(!mergeMode && (nje<njunctions && ejunction[nje]->end - currentstart < junctionsupport) &&
					(fs<0 || (uint)lstart[fs].predno>=ejunction[nje]->end) &&  // I do not want to miss any hard starts/ends
					(fe<0 || (uint)lend[fe].predno>=ejunction[nje]->end)) { // there is a junction ending soon here
				float covleft=get_cov(1,currentstart-refstart,ejunction[nje]->end-1-refstart,bpcov);
				float covright=get_cov(1,ejunction[nje]->end-refstart,2*ejunction[nje]->end - currentstart-1-refstart,bpcov);
				if(covleft<covright*(1-ERROR_PERC)) { // adjust start here if needed
					currentstart=ejunction[nje]->end;
					// I have to check ending junctions here again
					while(nje<njunctions && ejunction[nje]->end<=currentstart) { // read all junction ends at or before the current start -> assuming there are any (at this point, smaller junction ends should not be relevant to this bundle/currentstart
						if(ejunction[nje]->end==currentstart && (ejunction[nje]->strand+1) == 2*s) { // junction ends at current start and is on the same strand and not deleted
							end=1;
						}
						nje++;
					}
				}
			}
		}


		CGraphnode *graphnode=create_graphnode_unispg(s,g,currentstart,endbundle,graphno,bundlenode,bundle2graph,no2gnode); // creates a $graphno graphnode  with start at bundle start, and end at bundle end
		graphno++;

		// KH: For multi-samples, I should just skip this.
		/****************
		 **  original
		****************/
		if(end) { // I might have nodes finishing here; but I have a junction finishing here for sure
			//GStr cs((int)currentstart);
			//GVec<int> *e=ends[cs.chars()]; // HOW CAN I HAVE MORE THAN ONE NODE FINISHING HERE???; because this keeps all nodes that are linked by junctions here
			GVec<int> *e=ends[currentstart];
			if(e) {
				for(int i=0;i<e->Count();i++) {
					CGraphnode *node=no2gnode[s][g][e->Get(i)];
					node->child.Add(graphnode->nodeid);  // this node is the child of previous node
					graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
					// COUNT EDGE HERE
					edgeno++;
					fprintf(stderr,"1 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
				}
			}
			else { // I haven't seen nodes before that finish here (maybe due to error correction?) => link to source
				source->child.Add(graphnode->nodeid);  // this node is the child of source
				graphnode->parent.Add(source->nodeid); // this node has source as parent
				// COUNT EDGE HERE
				edgeno++;
				fprintf(stderr,"2 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
			}
		}
		else { // this node comes from source directly
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			// COUNT EDGE HERE
			edgeno++;
			fprintf(stderr,"3 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
		}
		/****************
		 **  original
		****************/

		bool completed=false;
		bool dropcov=false; // false(0) means start of bundle or junction end (raise in coverage); true(1) means junction start (drop in coverage)
		int nls=0; // index in longstart
		int nle=0; // index in longend

		do {
			while(nje<njunctions && (((int)ejunction[nje]->strand+1) != 2*s)) nje++; // skip junctions that don't have the same strand
			while(njs<njunctions && ((((int)junction[njs]->strand+1) != 2*s) || (junction[njs]->start<currentstart))) njs++; // junctions that start before the current graphnode and I haven't seen them before are part of a different bundle


			int minjunction = -1; // process next junction -> either a start or an ending whichever has the first position on the genome; if they have same position then process ending first
			if((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle))) {
				if(njs<njunctions && (junction[njs]->start<=endbundle) && junction[njs]->end>bundle_end) njs++;
				else {
					if(nje<njunctions) { // there are still junctions endings
						if(njs<njunctions) { // there are still junctions starting
							minjunction = junction[njs]->start >= ejunction[nje]->end ? 1 : 0; // one of them is clearly before the endbundle from the initial if
						}
						else minjunction = 1;
					}
					else minjunction = 0;
				}
			}

			// fprintf(stderr,"minjunction=%d\n",minjunction);
			if(nje<njunctions) fprintf(stderr,"Found junction:%d-%d(%d)\n",ejunction[nje]->start,ejunction[nje]->end,ejunction[nje]->strand);


			if(minjunction == 0 ) { // found a start junction here

				// add guide starts/ends first
				if(processguide) {
					while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
					if(nge<guideedge.Count()) {

						while(true) {

							while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

							if(nge>=guideedge.Count() || guideedge[nge].val>=junction[njs]->start) break;

							uint gstart=guideedge[nge].val;
							uint gend=junction[njs]->start;
							bool sourceguide=false;
							if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
							nge++;
							if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
							else if(guideedge[nge-1].endval<currentstart) continue;

							while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
							if(nge<guideedge.Count() && guideedge[nge].val<junction[njs]->start) gend=guideedge[nge].val;

							// I need to check there is no other trimming needed due to drops from longreads
							if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,gstart,nls,nle,dropcov,!sourceguide,lstart,lend,
									graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

							if(sourceguide)	{
								graphnode=source2guide(s,g,refstart,gstart,gend,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
								dropcov=false;
							}
							else {
								graphnode=guide2sink(s,g,refstart,gstart,gend,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
								dropcov=true;
							}

						}
					}
				}
				// if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes
				else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,junction[njs]->start,nls,nle,dropcov,true,lstart,lend,
							graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
				if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


				dropcov=true;
				// if no trimming required just set the end of the node
				graphnode->end=junction[njs]->start; // set the end of current graphnode to here; introduce smaller nodes if trimming is activated
				uint pos=junction[njs]->start;
				while(njs<njunctions && junction[njs]->start==pos ) { // remember ends here
					if((junction[njs]->strand+1) == 2*s) {
						//seenjunc++;
						if(mergeMode && (int)junction[njs]->end==refend) { // this node goes straight to sink
							sink->parent.Add(graphnode->nodeid); // graphnode is the parent of sink: check to see if I have a conflict with this
							edgeno++; // count edge here
						}
						else {
							//GStr je((int)junction[njs]->end);
							//GVec<int> *e=ends[je.chars()];
							GVec<int> *e=ends[junction[njs]->end];
							if(!e) {
								e = new GVec<int>();
								//ends.Add(je.chars(),e);
								ends.Add(junction[njs]->end, e);
							}
							e->Add(graphnode->nodeid);
						}
					}
					njs++;
				}

				if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
					// see if I should just skip node
					if(endbundle-pos<junctionsupport) {
						while(njs<njunctions && junction[njs]->strand+1 != 2*s) njs++;
						if(!mergeMode && (njs>=njunctions || junction[njs]->start > endbundle) && (nje>=njunctions || ejunction[nje]->end > endbundle)) { // there are no more junctions starting within this bundle
							float covleft=get_cov(1,2*pos-endbundle+1-refstart,pos-refstart,bpcov);
							float covright=get_cov(1,pos+1-refstart,endbundle-refstart,bpcov);
							if(covright<covleft*(1-ERROR_PERC)) { // adjust start here if needed
								completed=true;
							}
						}
					}

					if(!completed) {
						//fprintf(stderr,"create graph 2\n");
						// CGraphnode *nextnode = create_graphnode_unispg_cov(s,g,pos+1,endbundle,graphno,1,bundlenode,bundle2graph,no2gnode);
						CGraphnode *nextnode = create_graphnode_unispg(s,g,pos+1,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
						graphno++;
						graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
						nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode
						// COUNT EDGE HERE
						edgeno++;
						fprintf(stderr,"4 Edge %d-%d, edgeno=%d nextnode: %u-%u pos=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno,nextnode->start,nextnode->end,pos);
						graphnode=nextnode;
					}
				}
				else completed=true;
			}
			else if(minjunction == 1) { // found a junction end here

				uint pos=ejunction[nje]->end;
				while(nje<njunctions && ejunction[nje]->end==pos) { // read all junction ends at the current start
					nje++;
				}
				if(graphnode->start<pos) { // last created node starts before the position of the new node I want to create

					// add guide starts/ends first
					if(processguide) {
						while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
						if(nge<guideedge.Count()) {

							while(true) {

								while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

								if(nge>=guideedge.Count() || guideedge[nge].val>=pos-1) break;

								uint start=guideedge[nge].val;
								uint end=pos-1;
								bool sourceguide=false;
								if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
								nge++;
								if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
								else if(guideedge[nge-1].endval<currentstart) continue;

								while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
								if(nge<guideedge.Count() && guideedge[nge].val<pos-1) end=guideedge[nge].val;

								// I need to check there is no other trimming needed due to drops from longreads
								if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,start,nls,nle,dropcov,!sourceguide,lstart,lend,
										graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

								if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
								else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

							}
						}
					}
					//if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes
					else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,pos-1,nls,nle,dropcov,false,lstart,lend,
								graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
					if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


					graphnode->end=pos-1; // set end of current graphnode here
					dropcov=false;
					// fprintf(stderr,"create graph 3\n");
					// CGraphnode *nextnode = create_graphnode_unispg_cov(s,g,pos,endbundle,graphno,1,bundlenode,bundle2graph,no2gnode);
					CGraphnode *nextnode = create_graphnode_unispg(s,g,pos,endbundle,graphno,bundlenode,bundle2graph,no2gnode);
					graphno++;
					graphnode->child.Add(nextnode->nodeid); // make nextnode a child of current graphnode
					nextnode->parent.Add(graphnode->nodeid);// make graphnode a parent of nextnode

					// COUNT EDGE HERE
					edgeno++;
					fprintf(stderr,"5 Edge %d-%d, edgeno=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno);

					graphnode=nextnode;
				}
				// GStr spos((int)pos);
				// GVec<int> *e=ends[spos.chars()]; // WHY DOESN'T THIS REPEAT THE SAME THING IN CASE THE START HASN'T BEEN ADJUSTED? because nje is bigger now than the ones that end at the currentstart

				GVec<int> *e=ends[pos];
				if(e) for(int i=0;i<e->Count();i++) {
					CGraphnode *node=no2gnode[s][g][e->Get(i)];
					node->child.Add(graphnode->nodeid);  // this node is the child of previous node
					graphnode->parent.Add(node->nodeid); // this node has as parent the previous node
					// COUNT EDGE HERE
					edgeno++;
					fprintf(stderr,"6 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
				}
			}
		} while((nje<njunctions && (ejunction[nje]->end<=endbundle)) || (njs<njunctions && (junction[njs]->start<=endbundle)));

		if(!completed) { // I did not finish node --> this will be an ending node

			// add guide starts/ends first
			if(processguide) {
				while(nge<guideedge.Count() && guideedge[nge].val<=graphnode->start) nge++;
				if(nge<guideedge.Count()) {

					while(true) {

						while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;

						if(nge>=guideedge.Count() || guideedge[nge].val>=endbundle) break;

						uint start=guideedge[nge].val;
						uint end=endbundle;
						bool sourceguide=false;
						if(guideedge[nge].val<guideedge[nge].endval) sourceguide=true;
						nge++;
						if(sourceguide) { if(guideedge[nge-1].endval>endbundle) continue;}
						else if(guideedge[nge-1].endval<currentstart) continue;

						while(nge<guideedge.Count() && guideedge[nge].strand!=s) nge++;
						if(nge<guideedge.Count() && guideedge[nge].val<endbundle) end=guideedge[nge].val;

						// I need to check there is no other trimming needed due to drops from longreads
						if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,start,nls,nle,dropcov,!sourceguide,lstart,lend,
								graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);

						if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
						else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

					}
				}
			}
			// if(trim && !processguide && !mergeMode) graphnode=trimnode(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes
			else if(longreads && (lstart.Count() || lend.Count())) graphnode=longtrim(s,g,refstart,endbundle,nls,nle,dropcov,true,lstart,lend,
						graphnode,source,sink,futuretr,graphno,bpcov,bundlenode,bundle2graph,no2gnode,edgeno);
			if(trim && !longreads && !mergeMode) graphnode=trimnode_all(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes


			graphnode->end=endbundle;
			// COUNT EDGE HERE (this is an edge to sink)
			edgeno++;
			fprintf(stderr,"7 Edge to sink from %d, edgeno=%d\n",graphnode->nodeid,edgeno);
		}
	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)
    

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
			else edgeno++; // node has no children so it might get linked to sink in traverse_dfs function
		}
	}

	/*****************************
	 ** Step 5: 'prune_graph_nodes' function
	 ** 	here I know graphno => I can see if it's too big
	 *****************************/
	if(!mergeMode && graphno>allowed_nodes) { // TODO: define allowed_nodes as a default in stringtie.cpp that varies with the memory
		graphno=prune_graph_nodes(graphno,s,g,bundle2graph,bnode.Count(),no2gnode,junction,edgeno,futuretr,sink);
	}

	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

	if(mergeMode) { // I might have a bunch of sink's parents that are not linked to sink
		for(int i=0;i<sink->parent.Count();i++) {
			CGraphnode *node=no2gnode[s][g][sink->parent[i]];
			node->child.Add(sink->nodeid);
		}
	}

	// fprintf(stderr,"This graph has %d nodes and %d edges and starts at lastpos=%d\n",graphno,edgeno,graphno);
	lastgpos=graphno; // nodes are from 0 to graphno-1, so the first "available" position in GBitVec is graphno










	/***************************
	 ** In C++ visualization
 	 ***************************/
	// if(no2gnode[s][g].Count()) {
	// 	GStr strand_symbol;
	// 	if (s == 0) {
	// 		strand_symbol = "-";
	// 	} else if (s == 1) {
	// 		strand_symbol = "+";
	// 	}

	// 	fprintf(stderr,"Traversing the created graph!!!\n");
	// 	fprintf(stderr,"Digraph %d_%d_%d_%d {", bdata->start,bdata->end, s, g);
	// 	// graphno[s][b]: number of nodes in graph.
	// 	if(no2gnode[s][g].Count()) {

	// 		// for(int nd=0;nd<no2gnode[s][g].Count();nd++) {
	// 		for(int nd=1;nd<no2gnode[s][g].Count()-1;nd++) {
	// 			fprintf(stderr,"%d[start=%d end=%d cov=%f];",nd,no2gnode[s][g][nd]->start,no2gnode[s][g][nd]->end,no2gnode[s][g][nd]->cov);
	// 			exon_tmp.clear();
	// 			// exon_tmp.push_back(no2gnode[s][g][nd]->start);
	// 			// exon_tmp.push_back(no2gnode[s][g][nd]->end);
	// 			exon_tmp.push_back(nd*3);
	// 			exon_tmp.push_back(nd*3+1);
	// 			exonIntervals.push_back(exon_tmp);
	// 			// for (int i = no2gnode[s][g][nd]->start; i < no2gnode[s][g][nd]->end; i++) {
	// 			// 	fprintf(node_cov_bed, "chr22\t%d\t%d\tNODE\t%d\t%s\n", i, i+1, no2gnode[s][g][nd]->cov, strand_symbol.chars());
	// 			// }
	// 		}

	// 		for(int nd=0;nd<no2gnode[s][g].Count();nd++) {
	// 			// fprintf(stderr,"Node %d with parents:",i);
	// 			for(int c=0;c<no2gnode[s][g][nd]->child.Count();c++) {
	// 				fprintf(stderr,"%d->",nd);			
	// 				fprintf(stderr,"%d;",no2gnode[s][g][nd]->child[c]);
	// 				if (no2gnode[s][g][nd]->end == 0 || no2gnode[s][g][ no2gnode[s][g][nd]->child[c] ] -> start == 0) {

	// 				} else {
	// 					intron_tmp.clear();
	// 					// intron_tmp.push_back(no2gnode[s][g][nd]->end);
	// 					// intron_tmp.push_back(no2gnode[s][g][ no2gnode[s][g][nd]->child[c] ] -> start);
	// 					intron_tmp.push_back(nd*3+1);
	// 					intron_tmp.push_back(no2gnode[s][g][nd]->child[c]*3);
	// 					intronIntervals.push_back(intron_tmp);
	// 					// fprintf(edge_cov_bed, "chr22\t%d\t%d\tJUNCID\t%d\t%s\n", no2gnode[s][g][nd]->end, no2gnode[s][g][ no2gnode[s][g][nd]->child[c] ] -> start, 10, strand_symbol.chars());
	// 				}
	// 			}
	// 		}
	// 	}
	// 	fprintf(stderr,"}\n");
	// }

	// 	fprintf(stderr, "exonIntervals size: %d \n", exonIntervals.size());
	// 	fprintf(stderr, "intronIntervals size: %d \n", intronIntervals.size());
	// if (!exonIntervals.empty()) {
	// 	if (exonIntervals.size() >= 4) {
	// 		draw(exonIntervals, intronIntervals, string(plot_dir.chars())+"/plot_"+to_string(bundle_start)+"_"+to_string(bundle_end)+"_"+to_string(s)+"_"+to_string(g)+".png");
	// 	}
	// }
	/***************************
 	 ** In C++ visualization
 	 ***************************/







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

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"Add future transfrag[%d][%d]= %d with %d nodes n1=%d n2=%d graphno=%d, abundance=%f and pattern",s,g,transfrag[s][g].Count(),tr->nodes.Count(),n1,n2,graphno,futuretr[i+2]);
				//printBitVec(trpat);
				fprintf(stderr,"\n");
			}
			*/

			/*if(mixedMode) {
				tr->abundance*=2;
			}*/

			transfrag[s][g].Add(tr);
			if(mixedMode) {
				CTransfrag *longtr=new CTransfrag(nodes,trpat,futuretr[i+2]);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			}
			else
			if(longreads)// || mixedMode)
				transfrag[s][g].Last()->longread=true;
		}
	}


	/*****************************
	 ** Step 7: 'traverse_dfs_unispg' function
	 ** 	finished reading bundle -> now create the parents' and children's patterns
	 *****************************/
	GVec<bool> visit;
	visit.Resize(graphno);
	GBitVec parents(graphno+edgeno);

	// fprintf(stderr,"traverse graph[%d][%d] now with %d nodes, %d edges and lastgpos=%d....\n",s,g,graphno,edgeno,lastgpos);//edgeno=0;
	traverse_dfs_unispg(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
	// fprintf(stderr,"done traversing with edgeno=%d lastgpos=%d\n",edgeno,lastgpos);

	/*
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
	*/

	// delete variables created here, like e.g. ends; do I need to delete the GVec<int> elements created too?
	ends.Clear();

	return(graphno);
}




CTreePat *construct_treepat_unispg(int gno, GIntHash<int>& gpos,GPVec<CTransfrag>& transfrag) {

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













void process_transfrags_unispg(int s, int gno,int edgeno,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no,
		GIntHash<int> &gpos,GVec<CGuide>& guidetrf,GList<CPrediction>& pred,GVec<int>& trflong) {
	// /*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"\n***********************************************\n");
		fprintf(stderr,"*********** Inside process_transfrags_unispg: \n");
		fprintf(stderr,"***********************************************\n");
		fprintf(stderr,"There are %d transfrags before clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d](%f,%f) long=%d short=%d usepath=%d:",i,transfrag[i]->abundance,transfrag[i]->srabund,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	// */
	/****************************************
	 ** Step 1:
	 ** add all guide patterns to the set of transfrags so that I can have a "backbone" for each guide
	 ** I need this because there might be an incompatible transfrag connecting the nodes in the guide 
	 ****************************************/
	//fprintf(stderr,"There are %d guides\n",guidetrf.Count());
	for(int i=0;i<guidetrf.Count();i++) if(guidetrf[i].trf->guide){

		CTransfrag *t=NULL;
		bool add=true;
		if(longreads || mixedMode) {
			/*guidetrf[i].trf->pattern[0]=0;
			guidetrf[i].trf->pattern[gno-1]=0;
			int *pos=gpos[edge(0,guidetrf[i].trf->nodes[1],gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=0;
			pos=gpos[edge(guidetrf[i].trf->nodes[guidetrf[i].trf->nodes.Count()-2],guidetrf[i].trf->nodes.Last(),gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=0;
			guidetrf[i].trf->nodes.Pop();
			guidetrf[i].trf->nodes.Shift();*/

			t=findtrf_in_treepat(gno,gpos,guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,tr2no); // I need to adjust first/last node
			if(!t) { // t is NULL
				float abund=0;
				if(mixedMode) abund=trthr*ERROR_PERC;
				t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,abund);
				t->longread=true;
			}
			else add=false;
		}
		else {
			t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,trthr*ERROR_PERC);
		}

		if(!longreads) {
			if(includesource) {
				guidetrf[i].trf->nodes.Insert(0,0); // I need to comment this if I need path not to include the source
				guidetrf[i].trf->pattern[0]=1;
				int *pos=gpos[edge(0,guidetrf[i].trf->nodes[1],gno)];
				if(pos) guidetrf[i].trf->pattern[*pos]=1;
			}
			int sink=gno-1;
			guidetrf[i].trf->nodes.Add(sink);
			guidetrf[i].trf->pattern[sink]=1;
			int *pos=gpos[edge(guidetrf[i].trf->nodes[guidetrf[i].trf->nodes.Count()-2],guidetrf[i].trf->nodes.Last(),gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=1;
		}



		/*
		float abund=0;
		if(!longreads) {
			abund=trthr*ERROR_PERC;
		}
		CTransfrag *t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,abund);
		*/
		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Add guidetrf with nodes:");
			for(int j=0;j<guidetrf[i].trf->nodes.Count();j++) fprintf(stderr," %d",guidetrf[i].trf->nodes[j]);
			//fprintf(stderr," and pattern: ");
			//printBitVec(guidetrf[i].trf->pattern);
			fprintf(stderr,"\n");
			fprintf(stderr,"Found transfrag %p with nodes:",t);
			for(int j=0;j<t->nodes.Count();j++) fprintf(stderr," %d",t->nodes[j]);
			//fprintf(stderr," and pattern: ");
			//printBitVec(guidetrf[i].trf->pattern);
			fprintf(stderr,"\n");

		}
		*/

		/*if(!longreads) {
			t->pattern[0]=0;
			t->pattern[gno-1]=0;
			int *pos=gpos[edge(0,t->nodes[1],gno)];
			if(pos) t->pattern[*pos]=0;
			pos=gpos[edge(t->nodes[t->nodes.Count()-2],t->nodes.Last(),gno)];
			if(pos) t->pattern[*pos]=0;
			if(add) {
				t->nodes.Pop();
				t->nodes.Shift();
			}
		}*/
		t->guide=1+guidetrf[i].g;
		//fprintf(stderr,"t->guide set to=%d\n",1+guidetrf[i].g);
		t->longstart=no2gnode[t->nodes[0]]->start;
		t->longend=no2gnode[t->nodes.Last()]->end;
		//if(longreads) t->usepath=guidetrf[i].g; // guide index
		no2gnode[t->nodes[0]]->hardstart=1;  // I can always trust a guide's start
		no2gnode[t->nodes.Last()]->hardend=1; // I can always trust a guide's end
		if(add) transfrag.Add(t);
	}
	else if(!longreads) {
		if(includesource) {
			guidetrf[i].trf->nodes.Insert(0,0); // I need to comment this if I need path not to include the source
			guidetrf[i].trf->pattern[0]=1;
			int *pos=gpos[edge(0,guidetrf[i].trf->nodes[1],gno)];
			if(pos) guidetrf[i].trf->pattern[*pos]=1;
		}
		int sink=gno-1;
		guidetrf[i].trf->nodes.Add(sink);
		guidetrf[i].trf->pattern[sink]=1;
		int *pos=gpos[edge(guidetrf[i].trf->nodes[guidetrf[i].trf->nodes.Count()-2],guidetrf[i].trf->nodes.Last(),gno)];
		if(pos) guidetrf[i].trf->pattern[*pos]=1;
	}

	GPVec<CTransfrag> srfrag(false);
	/****************************************
	 ** Step 2:
	 ** eliminate transfrags below threshold (they represent noise) if they don't come from source
	 ****************************************/
	if(!eliminate_transfrags_under_thr(gno,gpos,transfrag,tr2no,trthr,srfrag) && srfrag.Count()) { // long transfrags but only to source/sink
		srfrag.Clear();
	}

	// /*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d] abund=%f long:%d: short:%d usepath:%d",i,transfrag[i]->abundance,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	// */

	GBitVec allpat(gno+edgeno);

	bool trsort=true;

	if(longreads) { // add source/sink links but only if they need to be added to explain the traversals in the graph
		transfrag.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete
		trsort=false;
		int source=0;
		int sink=gno-1;
		GVec<int> hassource(gno,-1); // remembers transcript number that links given node to source
		GVec<int> hassink(gno,-1); // remembers transcript number that links given node to sink
		GBitVec keepsource(gno); // if not set then I can remove link from node to source; keeps source link if it exists otherwise otherwise
		GBitVec keepsink(gno); // if not set then I can remove link from node to sink; keeps sink link if it exists otherwise otherwise
		GVec<CLongTrf> keeptrf; // keeps all potential transfrags that will be kept from most abundant to least, unassembled
		float zero=0;
		GVec<float> addsource(gno,zero);
		GVec<float> addsink(gno,zero);
		int edgedist=CHI_WIN; // I need to be consistent (if I change here then I need to change in update_abundance too)
		int ssdist=longintronanchor;

		/*
		{ // DEBUG ONLY
			//printTime(stderr);
			fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
			for(int i=0;i<transfrag.Count();i++) {
				fprintf(stderr,"transfrag[%d] abund=%f:",i,transfrag[i]->abundance);
				for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
				if(transfrag[i]->guide) fprintf(stderr," guide");
				fprintf(stderr,"\n");
			}
		}
		*/

		for(int t1=0;t1<transfrag.Count();t1++) {
			/*fprintf(stderr,"Consider t=%d with abund=%f guide=%d and nodes:",t1,transfrag[t1]->abundance,transfrag[t1]->guide);
			for(int j=0;j<transfrag[t1]->nodes.Count();j++) {
				if(j) {
					int *pos=gpos[edge(transfrag[t1]->nodes[j-1],transfrag[t1]->nodes[j],gno)];
					if(pos && transfrag[t1]->pattern[*pos]) {
						fprintf(stderr,"-");
					}
					else fprintf(stderr," ");
				}
				fprintf(stderr,"%d",transfrag[t1]->nodes[j]);
			} fprintf(stderr,"\n");*/
			if(!transfrag[t1]->nodes[0]) {
				hassource[transfrag[t1]->nodes[1]]=t1;
				//fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",transfrag[t1]->nodes[1],t1,transfrag[t1]->abundance);
			}
			else if(transfrag[t1]->nodes.Last()==gno-1) {
				hassink[transfrag[t1]->nodes[0]]=t1;
				//fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",transfrag[t1]->nodes[0],t1,transfrag[t1]->abundance);
			}
			else {
				if(eonly && !transfrag[t1]->guide) continue; // do not remember transfrags that are not guides
				if(!keepsource[transfrag[t1]->nodes[0]]) {
					if(transfrag[t1]->longstart) keepsource[transfrag[t1]->nodes[0]]=1; //fprintf(stderr,"keep source %d\n",transfrag[t1]->nodes[0]);}
					else if(no2gnode[transfrag[t1]->nodes[0]]->hardstart) keepsource[transfrag[t1]->nodes[0]]=1;
				}
				if(!keepsink[transfrag[t1]->nodes.Last()]) {
					if(transfrag[t1]->longend) keepsink[transfrag[t1]->nodes.Last()]=1;//fprintf(stderr,"keep sink %d\n",transfrag[t1]->nodes.Last());}
					else if(no2gnode[transfrag[t1]->nodes.Last()]) keepsink[transfrag[t1]->nodes.Last()]=1;
				}
				bool included=false;
				// a transfrag that starts at source and ends at sink can never be included in a kept transfrag, so I am safe to do next
				for(int t2=0; t2<keeptrf.Count();t2++) {
					int t[2]={t1,keeptrf[t2].t}; // t1 current, t2 the one I kept
					int len[4]={MAX_NODE,MAX_NODE,MAX_NODE,MAX_NODE};
					int ret=compatible_long(t,len,transfrag,no2gnode,gno,gpos);
					//fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
					if(ret){
						switch(ret) {
						case 1: // t[0] includes t[1]: it extends with introns on on one or both sides -> keep unless I can eliminate a previous one
							if(!transfrag[t[1]]->guide && transfrag[t1]->longstart && transfrag[t1]->longend && // t[1] might be included in t[0] so I might eliminate if it doesn't pass threshold
									(!no2gnode[transfrag[t[1]]->nodes[0]]->hardstart || transfrag[t[0]]->nodes[0] == transfrag[t[1]]->nodes[0] ) &&
									(!no2gnode[transfrag[t[1]]->nodes.Last()]->hardend || transfrag[t[0]]->nodes.Last() == transfrag[t[1]]->nodes.Last())) {
								//if(transfrag[t[0]]->abundance>(1-ERROR_PERC/DROP)*transfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
								if(transfrag[t[0]]->abundance>DROP*transfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
									if(len[1]<ssdist && len[3]<ssdist) { // prefer t[0] instead of t[1]
										keeptrf[t2].t=t1;
										keeptrf[t2].cov+=transfrag[t1]->abundance;
										keeptrf[t2].group.Add(t1);
										included=true; // I do not want to store transcript
										//fprintf(stderr,"trf %d includes %d\n",t[0],t[1]);
									}
								}
							}
							break;
						case 2: // t[1] includes t[0]: extends with introns past ends of t[0] (t[1] possibly includes t[0]); t[1] is more abundant than t0
							//if(transfrag[t[1]]->guide || transfrag[t[1]]->abundance>(1-ERROR_PERC/DROP)*transfrag[t[0]]->abundance) {
							if(!transfrag[t[0]]->guide &&
									(!no2gnode[transfrag[t[0]]->nodes[0]]->hardstart || transfrag[t[0]]->nodes[0] == transfrag[t[1]]->nodes[0]) &&
									(!no2gnode[transfrag[t[0]]->nodes.Last()]->hardend || transfrag[t[0]]->nodes.Last() == transfrag[t[1]]->nodes.Last())) {
								if(len[1]<ssdist && len[3]<ssdist) {
									keeptrf[t2].cov+=transfrag[t1]->abundance;
									keeptrf[t2].group.Add(t1);
									included=true;
									//fprintf(stderr,"trf %d is intronic including %d\n",t[1],t[0]);
								}
							}
							//}
							break;
						case 3: // t1 and t0 are compatible --> just look for the edges; t1 goes further apart
							if(transfrag[t[0]]->nodes[0]!=transfrag[t[1]]->nodes[0] && transfrag[t[0]]->nodes.Last()!=transfrag[t[1]]->nodes.Last() &&
									((no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && !no2gnode[transfrag[t[1]]->nodes[0]]->hardstart &&
											!no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && no2gnode[transfrag[t[1]]->nodes.Last()]->hardend) ||
											(!no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[1]]->nodes[0]]->hardstart &&
													no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && !no2gnode[transfrag[t[1]]->nodes.Last()]->hardend))) {
								// these two transcripts both have one good start and one different --> keep them both (different option would be to add them to another transfrag that is compatible and has both hardends but then it's more complicated
								break;
							}

							// I keep both if both are guides
							//if((!transfrag[t[0]]->guide || !transfrag[t[1]]->guide) && abs(len[0])<edgedist && abs(len[2])<edgedist) { // close by
							if((!transfrag[t[0]]->guide || !transfrag[t[1]]->guide) && len[0]<edgedist && len[2]<edgedist) { // close by
								if(transfrag[t[0]]->guide || (!transfrag[t[1]]->guide && no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[0]]->nodes.Last()]->hardend))
									keeptrf[t2].t=t1; // t[0] to replace t[1]
								keeptrf[t2].cov+=transfrag[t1]->abundance;
								keeptrf[t2].group.Add(t1);
								included=true;
								//fprintf(stderr,"trf %d %d equivalent start/ends\n",t[1],t[0]);
							}
							break;
						}
						if(included) break; // break from for loop
					}
				}

				if(!included){
					if(transfrag[t1]->guide || ((transfrag[t1]->longstart || no2gnode[transfrag[t1]->nodes[0]]->hardstart)&&
							(transfrag[t1]->longend || no2gnode[transfrag[t1]->nodes.Last()]->hardend))) { // if this is not included and has correct start/end
						CLongTrf kt(t1,transfrag[t1]->abundance);
						keeptrf.Add(kt);
						keeptrf.Last().group.Add(t1);
						//fprintf(stderr,"keep transfrag %d\n",t1);
					}
					else { // incomplete transcript, possibly wrong
						transfrag[t1]->weak=1;
						//fprintf(stderr,"Incomplete transcript %d\n",t1);
					}
				}
			}
		}


		char sign='-';
		if(s) { sign='+';}
		//GBitVec guidesource(gno);
		//GBitVec guidesink(gno);
		//for(int i=0;i<keeptrf.Count();i++) {
		for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier

			//fprintf(stderr,"Build source/sink for transfrag %d\n",keeptrf[i].t);
			int n1=transfrag[keeptrf[i].t]->nodes[0];
			int n2=transfrag[keeptrf[i].t]->nodes.Last();
			//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for transfrag %d\n",keeptrf[i].t);

			if(!rawreads && !transfrag[keeptrf[i].t]->guide && ((!no2gnode[n1]->hardstart && (hassource[n1]<0 || !keepsource[n1])) ||
					(!no2gnode[n2]->hardend && (hassink[n2]<0 || !keepsink[n2])))) {
			//if(!rawreads && (no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))) {
				trflong.Add(keeptrf[i].t);
			}

			/******* previous implementation here

			addsource[n1]=keeptrf[i].cov;
			addsink[n2]=keeptrf[i].cov;
			 *******/

			/*
			if(!addsource[n1] && hassource[n1]<0) {
				int startpos=no2gnode[n1]->start-refstart;
				if(startpos-CHI_THR<0 || startpos+CHI_THR>bpcov->Count()) addsource[n1]=1;
				else {
					addsource[n1]=(get_cov(1,startpos,startpos+CHI_THR-1,bpcov)-get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov)-
							get_cov(1,startpos-CHI_THR,startpos-1,bpcov)+get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov))/(DROP*CHI_THR);
				}
			}

			if(!addsink[n2]  && hassink[n1]<0) {
				int endpos=no2gnode[n2]->end-refstart;
				if(endpos-CHI_THR<0 || endpos+CHI_THR>bpcov->Count()-1) addsink[n2]=1;
				else {
					addsink[n2]=(get_cov(1,endpos-CHI_THR+1,endpos,bpcov)-get_cov(2-2*s,endpos-CHI_THR+1,endpos,bpcov)-
							get_cov(1,endpos+1,endpos+CHI_THR,bpcov)+get_cov(2-2*s,endpos+1,endpos+CHI_THR,bpcov))/(DROP*CHI_THR);
				}
			}*/

			// all
			for(int j=0;j<keeptrf[i].group.Count();j++) {
				if(n1==transfrag[keeptrf[i].group[j]]->nodes[0]) {
					addsource[n1]+=transfrag[keeptrf[i].group[j]]->abundance;
					//fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],transfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
				}
				if(n2==transfrag[keeptrf[i].group[j]]->nodes.Last()) {
					addsink[n2]+=transfrag[keeptrf[i].group[j]]->abundance;
					//fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],transfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
				}
			}


			if(rawreads) {
				GVec<GSeg> exons;
				int j=0;
				int len=0;
				int t=keeptrf[i].t;
				while(j<transfrag[t]->nodes.Count()) {
					int nodestart=no2gnode[transfrag[t]->nodes[j]]->start;
					int nodeend=no2gnode[transfrag[t]->nodes[j]]->end;
					len+=nodeend-nodestart+1;
					while(j+1<transfrag[t]->nodes.Count() && no2gnode[transfrag[t]->nodes[j]]->end+1==no2gnode[transfrag[t]->nodes[j+1]]->start) {
						j++;
						len+=no2gnode[transfrag[t]->nodes[j]]->len();
						nodeend=no2gnode[transfrag[t]->nodes[j]]->end;
					}
					GSeg exon(nodestart,nodeend);
					exons.Add(exon);
					j++;
				}
				uint tstart=exons[0].start;
				uint tend=exons.Last().end;
				if(transfrag[t]->longstart>tstart) {
					len-=transfrag[t]->longstart-tstart;
					tstart=transfrag[t]->longstart;
				}
				if(transfrag[t]->longend && transfrag[t]->longend<tend) {
					len-=tend-transfrag[t]->longend;
					tend=transfrag[t]->longend;
				}

				CPrediction *p=new CPrediction(s, NULL, tstart, tend, keeptrf[i].cov, sign, len);
				exons[0].start=tstart;
				exons.Last().end=tend;
				p->exons=exons;
				pred.Add(p);
			}
		}
		for(int i=keeptrf.Count()-1;i>=0;i--) {
			int n1=transfrag[keeptrf[i].t]->nodes[0];
			int n2=transfrag[keeptrf[i].t]->nodes.Last();

			if(transfrag[keeptrf[i].t]->guide || ((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))))
				trflong.Add(keeptrf[i].t);

		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"%d keeptrf:\n",keeptrf.Count());
			for(int i=0;i<keeptrf.Count();i++) {
				fprintf(stderr,"(%d) %d abund=%f keepcov=%f",i,keeptrf[i].t,transfrag[keeptrf[i].t]->abundance,keeptrf[i].cov);
				for(int j=0;j<transfrag[keeptrf[i].t]->nodes.Count();j++) {
					fprintf(stderr," %d",transfrag[keeptrf[i].t]->nodes[j]);
				}
				fprintf(stderr,"\n");

			}
			fprintf(stderr,"%d trflong:",trflong.Count());
			for(int i=0;i<trflong.Count();i++)
				fprintf(stderr," %d",trflong[i]);
			fprintf(stderr,"\n");
		}
		*/


		// add source/sink connections
		for(int i=1;i<gno-1;i++) {
			//fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f\n",i,hassource[i],addsource[i],hassink[i],addsink[i]);
			if(hassource[i]<0) {
				if(addsource[i] && no2gnode[i]->hardstart) {
					/*float abund=trthr;
					if(no2gnode[i]->hardstart) { // node i doesn't have source as parent but it should
						abund=addsource[i];
					}*/
					// else I am not that confident this is a start
					no2gnode[i]->parent.Insert(0,zero);
					no2gnode[0]->child.Add(i);
					GVec<int> nodes;
					nodes.Add(source);
					nodes.Add(i);
					//CTransfrag *t=new CTransfrag(nodes,allpat,abund);
					CTransfrag *t=new CTransfrag(nodes,allpat,addsource[i]);
					t->pattern[source]=1;
					t->pattern[i]=1;
					t->longread=true;
					transfrag.Add(t);
					//fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
				}
			}
			else {
				if(!keepsource[i]) { // it was a mistake to link this to source --> remove
					/*fprintf(stderr,"delete abundance of trf:");
					for(int j=0;j<transfrag[hassource[i]]->nodes.Count();j++) fprintf(stderr," %d",transfrag[hassource[i]]->nodes[j]);fprintf(stderr,"\n");*/
					transfrag[hassource[i]]->abundance=0;
				}
				else { // this one has source and it should be kept => update trf abundance to a more realistic value
					//fprintf(stderr,"source %d:%d abund=%f addsource=%f\n",i,no2gnode[i]->start,transfrag[hassource[i]]->abundance,addsource[i]);
					//if(guidesource[i] && transfrag[hassource[i]]->abundance<addsource[i])
					transfrag[hassource[i]]->abundance=addsource[i];
				}
			}

			if(hassink[i]<0) {
				if(addsink[i] && no2gnode[i]->hardend) {
					/*float abund=trthr;
					if(no2gnode[i]->hardend) {
						abund=addsink[i];
					}*/
					// else  not very confident about this end

					no2gnode[i]->child.Add(sink);
					no2gnode[sink]->parent.Add(i);
					GVec<int> nodes;
					nodes.Add(i);
					nodes.Add(sink);
					CTransfrag *t=new CTransfrag(nodes,allpat,addsink[i]);
					t->pattern[sink]=1;
					t->pattern[i]=1;
					t->longread=true;
					transfrag.Add(t);
					//fprintf(stderr,"Add link to sink: %d-%d with abundance=%f\n",i,sink,addsink[i]);
				}
			}
			else {
				if(!keepsink[i]) { // it was a mistake to link this to sink --> remove
					/*fprintf(stderr,"delete abundance of trf:");
					for(int j=0;j<transfrag[hassink[i]]->nodes.Count();j++) fprintf(stderr," %d",transfrag[hassink[i]]->nodes[j]);fprintf(stderr,"\n");*/
					transfrag[hassink[i]]->abundance=0;
				}
				else {
					//fprintf(stderr,"sink %d:%d abund=%f addsink=%f\n",i,no2gnode[i]->end,transfrag[hassink[i]]->abundance,addsink[i]);
					//if(guidesink[i] && transfrag[hassink[i]]->abundance<addsink[i])
					transfrag[hassink[i]]->abundance=addsink[i];
				}
			}
		}
	}
	else if(srfrag.Count() && mixedMode) { // add source/sink links but only if they need to be added to explain the traversals in the graph
			srfrag.Sort(longtrCmp); // most abundant transfrags in the graph come first, then the ones with most nodes, then the ones more complete
			int source=0;
			int sink=gno-1;
			GVec<int> hassource(gno,-1); // remembers transcript number that links given node to source
			GVec<int> hassink(gno,-1); // remembers transcript number that links given node to sink
			GBitVec keepsource(gno); // if not set then I can remove link from node to source; keeps source link if it exists otherwise otherwise
			GBitVec keepsink(gno); // if not set then I can remove link from node to sink; keeps sink link if it exists otherwise otherwise
			GVec<CLongTrf> keeptrf; // keeps all potential transfrags that will be kept from most abundant to least, unassembled
			float zero=0;
			GVec<float> addsource(gno,zero);
			GVec<float> addsink(gno,zero);
			int edgedist=CHI_WIN; // I need to be consistent (if I change here then I need to change in update_abundance too)
			int ssdist=longintronanchor;
			int ntrf=0; // number of trflong transcripts

			/*
			{ // DEBUG ONLY
				//printTime(stderr);
				fprintf(stderr,"\nThere are %d longsrtransfrags after clean up:\n",srfrag.Count());
				for(int i=0;i<srfrag.Count();i++) {
					fprintf(stderr,"srfrag[%d] longstart=%d longend=%d abund=%f t=%d:",i,srfrag[i]->longstart,srfrag[i]->longend,srfrag[i]->abundance,int(srfrag[i]->usepath));
					for(int j=0;j<srfrag[i]->nodes.Count();j++) fprintf(stderr," %d",srfrag[i]->nodes[j]);
					fprintf(stderr,"\n");
				}
			}
			*/

			for(int t1=0;t1<srfrag.Count();t1++) {
				/*fprintf(stderr,"Consider t=%d with abund=%f and nodes:",t1,srfrag[t1]->abundance);
				for(int j=0;j<srfrag[t1]->nodes.Count();j++) {
					if(j) {
						int *pos=gpos[edge(srfrag[t1]->nodes[j-1],srfrag[t1]->nodes[j],gno)];
						if(pos && srfrag[t1]->pattern[*pos]) {
							fprintf(stderr,"-");
						}
						else fprintf(stderr," ");
					}
					fprintf(stderr,"%d",srfrag[t1]->nodes[j]);
				} fprintf(stderr,"\n");*/
				if(!srfrag[t1]->nodes[0]) {
					hassource[srfrag[t1]->nodes[1]]=t1;
					//fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",srfrag[t1]->nodes[1],t1,srfrag[t1]->abundance);
				}
				else if(srfrag[t1]->nodes.Last()==gno-1) {
					hassink[srfrag[t1]->nodes[0]]=t1;
					//fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",srfrag[t1]->nodes[0],t1,srfrag[t1]->abundance);
				}
				else {
					if(eonly && !srfrag[t1]->guide) continue; // do not remember transfrags that are not guides
					if(!keepsource[srfrag[t1]->nodes[0]]) {
						if(srfrag[t1]->longstart) keepsource[srfrag[t1]->nodes[0]]=1;
						else if(no2gnode[srfrag[t1]->nodes[0]]->hardstart) keepsource[srfrag[t1]->nodes[0]]=1;
					}
					if(!keepsink[srfrag[t1]->nodes.Last()]) {
						if(srfrag[t1]->longend) keepsink[srfrag[t1]->nodes.Last()]=1;
						else if(no2gnode[srfrag[t1]->nodes.Last()]) keepsink[srfrag[t1]->nodes.Last()]=1;
					}
					bool included=false;
					// a transfrag that starts at source and ends at sink can never be included in a kept transfrag, so I am safe to do next
					for(int t2=0; t2<keeptrf.Count();t2++) {
						int t[2]={t1,keeptrf[t2].t}; // t1 current, t2 the one I kept
						int len[4]={MAX_NODE,MAX_NODE,MAX_NODE,MAX_NODE};
						int ret=compatible_long(t,len,srfrag,no2gnode,gno,gpos);
						//fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
						if(ret){
							switch(ret) {
							case 1: // t[0] includes t[1]: it extends with introns on on one or both sides -> keep unless I can eliminate a previous one
								if(!srfrag[t[1]]->guide && srfrag[t1]->longstart && srfrag[t1]->longend && // t[1] might be included in t[0] so I might eliminate if it doesn't pass threshold
										(!no2gnode[srfrag[t[1]]->nodes[0]]->hardstart || srfrag[t[0]]->nodes[0] == srfrag[t[1]]->nodes[0] ) &&
										(!no2gnode[srfrag[t[1]]->nodes.Last()]->hardend || srfrag[t[0]]->nodes.Last() == srfrag[t[1]]->nodes.Last())) {
									//if(srfrag[t[0]]->abundance>(1-ERROR_PERC/DROP)*srfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
									if(srfrag[t[0]]->abundance>DROP*srfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
										if(len[1]<ssdist && len[3]<ssdist) { // prefer t[0] instead of t[1]
											keeptrf[t2].t=t1;
											keeptrf[t2].cov+=srfrag[t1]->abundance;
											keeptrf[t2].group.Add(t1);
											included=true; // I do not want to store transcript
											//fprintf(stderr,"trf %d includes %d\n",t[0],t[1]);
										}
									}
								}
								break;
							case 2: // t[1] includes t[0]: extends with introns past ends of t[0] (t[1] possibly includes t[0]); t[1] is more abundant than t0
								//if(srfrag[t[1]]->guide || srfrag[t[1]]->abundance>(1-ERROR_PERC/DROP)*srfrag[t[0]]->abundance) {
								if(!srfrag[t[0]]->guide &&
										(!no2gnode[srfrag[t[0]]->nodes[0]]->hardstart || srfrag[t[0]]->nodes[0] == srfrag[t[1]]->nodes[0]) &&
										(!no2gnode[srfrag[t[0]]->nodes.Last()]->hardend || srfrag[t[0]]->nodes.Last() == srfrag[t[1]]->nodes.Last())) {
									if(len[1]<ssdist && len[3]<ssdist) {
										keeptrf[t2].cov+=srfrag[t1]->abundance;
										keeptrf[t2].group.Add(t1);
										included=true;
										//fprintf(stderr,"trf %d is intronic including %d\n",t[1],t[0]);
									}
								}
								//}
								break;
							case 3: // t1 and t0 are compatible --> just look for the edges; t1 goes further apart
								if(srfrag[t[0]]->nodes[0]!=srfrag[t[1]]->nodes[0] && srfrag[t[0]]->nodes.Last()!=srfrag[t[1]]->nodes.Last() &&
										((no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && !no2gnode[srfrag[t[1]]->nodes[0]]->hardstart &&
												!no2gnode[srfrag[t[0]]->nodes.Last()]->hardend && no2gnode[srfrag[t[1]]->nodes.Last()]->hardend) ||
												(!no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && no2gnode[srfrag[t[1]]->nodes[0]]->hardstart &&
														no2gnode[srfrag[t[0]]->nodes.Last()]->hardend && !no2gnode[srfrag[t[1]]->nodes.Last()]->hardend))) {
									// these two transcripts both have one good start and one different --> keep them both (different option would be to add them to another srfrag that is compatible and has both hardends but then it's more complicated
									break;
								}

								// I keep both if both are guides
								//if((!srfrag[t[0]]->guide || !srfrag[t[1]]->guide) && abs(len[0])<edgedist && abs(len[2])<edgedist) { // close by
								if((!srfrag[t[0]]->guide || !srfrag[t[1]]->guide) && len[0]<edgedist && len[2]<edgedist) { // close by
									if(srfrag[t[0]]->guide || (!srfrag[t[1]]->guide && no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && no2gnode[srfrag[t[0]]->nodes.Last()]->hardend))
										keeptrf[t2].t=t1; // t[0] to replace t[1]
									keeptrf[t2].cov+=srfrag[t1]->abundance;
									keeptrf[t2].group.Add(t1);
									included=true;
									//fprintf(stderr,"trf %d %d equivalent start/ends\n",t[1],t[0]);
								}
								break;
							}
							if(included) break; // break from for loop
						}
					}

					if(!included){
						if(srfrag[t1]->guide || ((srfrag[t1]->longstart || no2gnode[srfrag[t1]->nodes[0]]->hardstart)&&
								(srfrag[t1]->longend || no2gnode[srfrag[t1]->nodes.Last()]->hardend))) { // if this is not included and has correct start/end
							CLongTrf kt(t1,srfrag[t1]->abundance);
							keeptrf.Add(kt);
							keeptrf.Last().group.Add(t1);
							//fprintf(stderr,"keep srfrag %d\n",t1);
						}
						else { // incomplete transcript, possibly wrong
							srfrag[t1]->weak=1;
							//fprintf(stderr,"Incomplete transcript %d\n",t1);
						}
					}
				}
			}


			//GBitVec guidesource(gno);
			//GBitVec guidesink(gno);
			//for(int i=0;i<keeptrf.Count();i++) {
			for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier

				//fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);
				int n1=srfrag[keeptrf[i].t]->nodes[0];
				int n2=srfrag[keeptrf[i].t]->nodes.Last();
				//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);

				if(!srfrag[keeptrf[i].t]->guide && ((!no2gnode[n1]->hardstart && (hassource[n1]<0 || !keepsource[n1])) ||
						(!no2gnode[n2]->hardend && (hassink[n2]<0 || !keepsink[n2])))) {
				//if((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))) {
					srfrag[keeptrf[i].t]->usepath=-2-ntrf; // order of transfrag
					//fprintf(stderr,"keeptrf[%d].t=%d srfrag.usepath=%d\n",i,keeptrf[i].t,t);
					ntrf++;
					//trflong.Add(keeptrf[i].t);
				}

				/******* previous implementation here

				addsource[n1]=keeptrf[i].cov;
				addsink[n2]=keeptrf[i].cov;
				 *******/

				/*
				if(!addsource[n1] && hassource[n1]<0) {
					int startpos=no2gnode[n1]->start-refstart;
					if(startpos-CHI_THR<0 || startpos+CHI_THR>bpcov->Count()) addsource[n1]=1;
					else {
						addsource[n1]=(get_cov(1,startpos,startpos+CHI_THR-1,bpcov)-get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov)-
								get_cov(1,startpos-CHI_THR,startpos-1,bpcov)+get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov))/(DROP*CHI_THR);
					}
				}

				if(!addsink[n2]  && hassink[n1]<0) {
					int endpos=no2gnode[n2]->end-refstart;
					if(endpos-CHI_THR<0 || endpos+CHI_THR>bpcov->Count()-1) addsink[n2]=1;
					else {
						addsink[n2]=(get_cov(1,endpos-CHI_THR+1,endpos,bpcov)-get_cov(2-2*s,endpos-CHI_THR+1,endpos,bpcov)-
								get_cov(1,endpos+1,endpos+CHI_THR,bpcov)+get_cov(2-2*s,endpos+1,endpos+CHI_THR,bpcov))/(DROP*CHI_THR);
					}
				}*/

				// all
				for(int j=0;j<keeptrf[i].group.Count();j++) {
					if(n1==srfrag[keeptrf[i].group[j]]->nodes[0]) {
						addsource[n1]+=srfrag[keeptrf[i].group[j]]->abundance;
						//fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
					}
					if(n2==srfrag[keeptrf[i].group[j]]->nodes.Last()) {
						addsink[n2]+=srfrag[keeptrf[i].group[j]]->abundance;
						//fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
						}
					}

			}
			for(int i=keeptrf.Count()-1;i>=0;i--) {
				int n1=srfrag[keeptrf[i].t]->nodes[0];
				int n2=srfrag[keeptrf[i].t]->nodes.Last();

				if(srfrag[keeptrf[i].t]->guide || ((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) &&
						(no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2])))) {
					//trflong.Add(keeptrf[i].t);
					srfrag[keeptrf[i].t]->usepath=-2-ntrf; // mark transfrag that it needs to be part of trflong
					ntrf++;
				}
			}

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"%d keeptrf:\n",keeptrf.Count());
				for(int i=0;i<keeptrf.Count();i++) {
					fprintf(stderr,"(%d) %d abund=%f keepcov=%f",i,keeptrf[i].t,srfrag[keeptrf[i].t]->abundance,keeptrf[i].cov);
					for(int j=0;j<srfrag[keeptrf[i].t]->nodes.Count();j++) {
						fprintf(stderr," %d",srfrag[keeptrf[i].t]->nodes[j]);
					}
					fprintf(stderr,"\n");

				}
				fprintf(stderr,"%d trflong:",trflong.Count());
				for(int i=0;i<trflong.Count();i++)
					fprintf(stderr," %d",trflong[i]);
				fprintf(stderr,"\n");
			}
			*/


			// add source/sink connections
			for(int i=1;i<gno-1;i++) {
				//fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f hardstart=%d hardend=%d\n",i,hassource[i],addsource[i],hassink[i],addsink[i],no2gnode[i]->hardstart,no2gnode[i]->hardend);
				if(hassource[i]<0) {
					if(addsource[i] && no2gnode[i]->hardstart) {
						/*float abund=trthr;
						if(no2gnode[i]->hardstart) { // node i doesn't have source as parent but it should
							abund=addsource[i];
						}*/
						// else I am not that confident this is a start
						no2gnode[i]->parent.Insert(0,zero);
						no2gnode[0]->child.Add(i);
						GVec<int> nodes;
						nodes.Add(source);
						nodes.Add(i);
						//CTransfrag *t=new CTransfrag(nodes,allpat,abund);
						CTransfrag *t=new CTransfrag(nodes,allpat,addsource[i]);
						t->pattern[source]=1;
						t->pattern[i]=1;
						t->longread=true; /// this is only true for longreads;
						transfrag.Add(t);
						//fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
					}
				}
				else {
					/*if(keepsource[i]) {
						srfrag[hassource[i]]->abundance+=addsource[i];
					}*/
					if(!keepsource[i]) { // it was a mistake to link this to source --> remove
						/*fprintf(stderr,"delete abundance of trf:");
						for(int j=0;j<srfrag[hassource[i]]->nodes.Count();j++) fprintf(stderr," %d",srfrag[hassource[i]]->nodes[j]);fprintf(stderr,"\n");*/
						srfrag[hassource[i]]->abundance=0;
					}
					else { // this one has source and it should be kept => update trf abundance to a more realistic value
						//fprintf(stderr,"source %d:%d abund=%f addsource=%f\n",i,no2gnode[i]->start,srfrag[hassource[i]]->abundance,addsource[i]);
						//if(guidesource[i] && srfrag[hassource[i]]->abundance<addsource[i])
						srfrag[hassource[i]]->abundance=addsource[i];
					}
				}

				if(hassink[i]<0) {
					if(addsink[i] && no2gnode[i]->hardend) {
						/*float abund=trthr;
						if(no2gnode[i]->hardend) {
							abund=addsink[i];
						}*/
						// else  not very confident about this end

						no2gnode[i]->child.Add(sink);
						no2gnode[sink]->parent.Add(i);
						GVec<int> nodes;
						nodes.Add(i);
						nodes.Add(sink);
						CTransfrag *t=new CTransfrag(nodes,allpat,addsink[i]);
						t->pattern[sink]=1;
						t->pattern[i]=1;
						t->longread=true; //// this is only true for longreads
						transfrag.Add(t);
						//fprintf(stderr,"Add link to sink: %d-%d with abundance=%f\n",i,sink,addsink[i]);
					}
				}
				else {
					if(!keepsink[i]) { // it was a mistake to link this to sink --> remove
						/*fprintf(stderr,"delete abundance of trf:");
						for(int j=0;j<srfrag[hassink[i]]->nodes.Count();j++) fprintf(stderr," %d",srfrag[hassink[i]]->nodes[j]);fprintf(stderr,"\n");*/
						srfrag[hassink[i]]->abundance=0;
					}
					else {
						//fprintf(stderr,"sink %d:%d abund=%f addsink=%f\n",i,no2gnode[i]->end,srfrag[hassink[i]]->abundance,addsink[i]);
						//if(guidesink[i] && srfrag[hassink[i]]->abundance<addsink[i])
						srfrag[hassink[i]]->abundance=addsink[i];
					}
					/*if(keepsink[i]) {
						//fprintf(stderr,"Add %f to longtransfrag=%d\n",addsink[i],hassink[i]);
						srfrag[hassink[i]]->abundance+=addsink[i];
					}*/
				}
			}
			//fprintf(stderr,"There should be %d trflong transcripts\n",ntrf);
			trflong.Resize(ntrf,-1);
	}


	/****************************************
	 ** Step 3:
	 ** add edges between disconnected parent-child nodes
	 ****************************************/
	for(int t=0;t<transfrag.Count();t++) allpat=allpat | transfrag[t]->pattern;

	/****************************************
	 ** Step 4:
	 ** for all nodes check if there is a connection to child
	 ****************************************/
	for(int i=1;i<gno-1;i++) { 
		CGraphnode *n=no2gnode[i];
		for(int c=0;c<n->child.Count();c++) {
			int *pos=gpos[edge(i,n->child[c],gno)];
			if(pos && !allpat[*pos]) {
				GVec<int> nodes;
				nodes.Add(i);
				nodes.Add(n->child[c]);
				GBitVec trpat(gno+edgeno);
				trpat[i]=1;
				trpat[n->child[c]]=1;
				trpat[*pos]=1;
				CTransfrag *t=new CTransfrag(nodes,trpat,trthr);
				if(longreads) t->longread=true;
				transfrag.Add(t);
			}
		}
	}

	/****************************************
	 ** Step 5:
	 ** sort transfrag with smallest being the one that has the most nodes, and ties are decided by the abundance (largest abundance first); last transfrags all have 1 node
	 ****************************************/
	if(trsort)
		transfrag.Sort(trCmp);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags that remained\n",transfrag.Count());
	}
	*/

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after sorting:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d]:",i);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr," abundance=%f usepath=%.1f\n",transfrag[i]->abundance,transfrag[i]->usepath);
		}
	}
	*/


	GVec<int> incompletetrf; //remembers incomplete transfrags (the ones that don't have edges between two consecutive nodes


	/****************************************
	 ** Step 6:
	 ** create compatibilities
	 ****************************************/
	for(int t1=0;t1<transfrag.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations

		// update nodes
		int n1=transfrag[t1]->nodes.Count();

		if(mixedMode) {
			if(transfrag[t1]->usepath<-1) { // this is a mixedMode long transfrag (TODO: consider using this for longreads mode also)

				//fprintf(stderr,"Add transcript %d to trflong with value %d\n",t1,(int)transfrag[t1]->usepath);
				trflong[abs(2+(int)transfrag[t1]->usepath)]=t1;
			}
			//transfrag[t1]->usepath=-1; // restore order
			transfrag[t1]->usepath=transfrag[t1]->abundance; // store abundance for later
		}

		if(n1>1) { // add transfrag to nodes' in and out; if a transfrag only has one node then it is not added to a node; I might want to change this for the computation of fpkm's
			bool incomplete = false;
			bool nosplice=true; // try to give less priority to unspliced reads vs spliced reads
			for(int n=0;n<n1;n++) { // for all nodes in transfrag

				if(nosplice && n) { // reduce abundance of continuous transfrags
					if(transfrag[t1]->nodes[n]!=1+transfrag[t1]->nodes[n-1] || no2gnode[transfrag[t1]->nodes[n]]->start-1!=no2gnode[transfrag[t1]->nodes[n-1]]->end) {
						nosplice=false;
					}
				}

				if(n && n<transfrag[t1]->nodes.Count()-1) {// not first or last node
					// add t1 to in and out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					if(transfrag[t1]->nodes[n-1] && transfrag[t1]->nodes[n]<gno-1) {
						// check if transfrag t1 is incomplete between node[n-1] and node [n]
						int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
						if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
							incomplete = assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
					}
				}
				else if(n) { // last but not first node
					// add t1 to in of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					if(transfrag[t1]->nodes[n-1] && transfrag[t1]->nodes[n]<gno-1) {
						// check if transfrag t1 is incomplete between node[n-1] and node [n]
						int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];
						if(!pos || !transfrag[t1]->pattern[*pos]) // there is no edge between node[n-1] and node[n]
							incomplete = assign_incomplete_trf_to_nodes(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
					}
				}
				else { // first node -> only add transfrag to out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);
				}
			}

			//if(nosplice) transfrag[t1]->abundance*=(1-isofrac);
			//transfrag[t1]->abundance*=0.5;


			if(incomplete) incompletetrf.Add(t1);
			else transfrag[t1]->real=true;

		}
		//else if(longreads) no2gnode[n1]->trf.Add(t1);
		/*
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}
		*/
	} // end for(int t1=0;t1<transfrag.Count();t1++)


	if(srfrag.Count() && !mixedMode) {
		srfrag.Sort(trCmp); // always start with largest super-read to solve
		for(int u=0;u<srfrag.Count();u++) //process_srfrag(srfrag[u],transfrag,no2gnode,gno,gpos);
		  if(!srfrag[u]->abundance) srfrag[u]->abundance=srfrag[u]->srabund*ERROR_PERC;
	}


	/****************************************
	 ** Step 7:
	 **  set source-to-child transfrag abundances: optional in order not to keep these abundances too low:
	 **  update the abundances of the transfrags coming in from source and going to a node that doesn't have other parents than source
	 **  * this part was removed to improve performance
	 ****************************************/
	CGraphnode *source=no2gnode[0];
	for(int i=0;i<source->child.Count();i++) {
		float abundance=0;
		int t0=-1;
		if(no2gnode[source->child[i]]->parent.Count()==1 && !no2gnode[source->child[i]]->parent[0]) { // source is the only parent of node
			for(int j=0;j<no2gnode[source->child[i]]->trf.Count();j++) {
				int t=no2gnode[source->child[i]]->trf[j];
				if(transfrag[t]->nodes.Last()==source->child[i]) t0=t;
				else abundance+=transfrag[t]->abundance;
			}
			if(t0>-1 && transfrag[t0]->abundance) { // found transfrag from source to node and the transfrag wasn't deleted
				transfrag[t0]->abundance=abundance;
			}
		}
	}
	// */

	for(int t=0;t<incompletetrf.Count();t++)
		transfrag[incompletetrf[t]]->real=trf_real(incompletetrf[t],no2gnode,transfrag,gpos,gno);

}



















int build_graphs_unispg(BundleData* bdata, int fidx) {
	// Clear the unispg_gp.
	// unispg_gp->Clear();
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
	 ** Step 0: this part is for setting guides for introns are covered by at least one read. 
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

				//fprintf(stderr,"Look to add guide g=%d start %d-%d and end %d-%d on strand %d\n",g,guides[g]->start,guides[g]->exons[0]->end,guides[g]->end,guides[g]->exons.Last()->start,s);

				int uses=s;
				if(s<0) uses=0;
				// guide edge
				GEdge ge(guides[g]->start,guides[g]->exons[0]->end,uses);
				int idx=guideedge.IndexOf(ge);

				//fprintf(stderr,"look for ge(%d,%d,%d) => start idx=%d\n",ge.val,ge.endval,ge.strand,idx);

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

				//fprintf(stderr,"look for ge(%d,%d,%d) => end idx=%d\n",ge.val,ge.endval,ge.strand,idx);

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

		/*
		fprintf(stderr,"check read[%d]:%d-%d:%d refstart=%d w/exons:",n,rd.start,rd.end,rd.strand,refstart);
		for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
		fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);
		i=0;
		*/

		while(i<rd.juncs.Count()) {
			CJunction& jd=*(rd.juncs[i]);
			// fprintf(stderr, " read junc %d-%d:%d nreads=%f nreads_good=%f nm=%f support=%f,%f between exons:%d-%d and %d-%d\n", jd.start, jd.end, jd.strand,jd.nreads,jd.nreads_good,jd.nm,jd.leftsupport,jd.rightsupport,rd.segs[i].start,rd.segs[i].end,rd.segs[i+1].start,rd.segs[i+1].end);

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


	/*****************************
	 ** Step 4: 'merge_fwd_groups' function
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
	 ** Step 5: form bundles here
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
	 ** Step 6: 'set_strandcol' function
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
	 ** Step 7: 'create_bundle' function. 
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
	 ** Step 8: Clean up no longer needed variables
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
	 ** Step 9: 'get_covered' function
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
	 ** Step 10: 'CPrediction': constructor
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
	 ** Step 11: Defining parameters here!!!!
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
	 	 **     'create_graph_unispg', 'construct_treepat_unispg'
		 *****************************/
    	for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions

        	// guides appear to be sorted by start --> CHECK THIS!!
        	int g=0;
        	int ng=guides.Count();

    		int s=sno/2; // adjusted strand due to ignoring neutral strand

			int graph_no = bundle[sno].Count();
			fprintf(stderr, "** graph_no: %d\n", graph_no);
			// bundle[sno].Count();

    		char strnd='-';
    		if(s) strnd='+';

    		bundle2graph[s]=NULL;
    		if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
    		transfrag[s]=NULL;
    		no2gnode[s] = NULL;
    		tr2no[s]=NULL;
    		gpos[s]=NULL;

			// fprintf(stderr, "0 bundle[sno].Count(): %d\n", bundle[sno].Count());
    		if(graph_no) {
    			transfrag[s]=new GPVec<CTransfrag>[graph_no]; // for each bundle I have a graph ? only if I don't ignore the short bundles

    			no2gnode[s]=new GPVec<CGraphnode>[graph_no];
    			gpos[s]=new GIntHash<int>[graph_no];

    			GCALLOC(tr2no[s],graph_no*sizeof(CTreePat *));
    			bno[s]=graph_no;

    			for(int b=0;b<graph_no;b++) {

    				graphno[s].cAdd(0);
    				edgeno[s].cAdd(0);
    				lastgpos[s].cAdd(0);
    				// I am overestmating the edgeno below, hopefully not by too much

    				//fprintf(stderr,"Bundle is: %d - %d start at g=%d sno=%d b=%d\n",bnode[sno][bundle[sno][b]->startnode]->start,bnode[sno][bundle[sno][b]->lastnodeid]->end,g,sno,b);

					fprintf(stderr, "bnode[%d][bundle[%d][%d]->startnode]->start: %d \n", sno, sno, b, bnode[sno][bundle[sno][b]->startnode]->start);

					fprintf(stderr, "bnode[%d][bundle[%d][%d]->lastnodeid]->end: %d \n", sno, sno, b, bnode[sno][bundle[sno][b]->lastnodeid]->end);

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
    						edgeno[s][b]+=2; // this is an overestimate: possibly I have both an extra source and an extra sink link
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
    					graphno[s][b]=create_graph_unispg(s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,trims); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure
    					*
    					*/

						/****************
						 **  KH modify
						 ****************/
    					// create graph then
    					graphno[s][b]=create_graph_unispg(refstart,s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,gpos,bdata,edgeno[s][b],lastgpos[s][b],guideedge); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

    					if(graphno[s][b]) tr2no[s][b]=construct_treepat_unispg(graphno[s][b],gpos[s][b],transfrag[s][b]);
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
				fprintf(stderr, ">>>>>>> 1 bundle[%d].Count(): %d\n", s, bno[s]);
				fprintf(stderr, ">>>>>>> 1 graphno[%d][%d]: %d\n", s, b, graphno[s][b]);
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

		// unispg_gp->WriteGraphGp();

		/*****************************
		 ** 3-1. Write out global splice graph in DOT format
		*****************************/
		/****************
		 **  KH Adding 
		 ****************/
		// if (universal_splice_graph) {
		// 	//  DOT file outut here 
		// 	//  not capacity and rate 
		// 	//  only edge weight
		// 	for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
		// 		int s=sno/2; // adjusted strand due to ignoring neutral strand
		// 		int g_idx = 0;
		// 		for(int b=0;b<bundle[sno].Count();b++) {
		// 			fprintf(stderr, "New writing out place!! Start writing out DOT file!!\n");
		// 			fprintf(stderr,"after traverse:\n");
		// 			if(graphno[s][b]) {
		// 				fprintf(uinigraph_out,"strict digraph %d_%d_%d_%d {", refstart, refend, s, g_idx);
		// 				// graphno[s][b]: number of nodes in graph.
		// 				if(graphno[s][b]) {
		// 					for(int nd=1;nd<graphno[s][b]-1;nd++)
		// 						fprintf(uinigraph_out,"%d[start=%d end=%d cov=%f];",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end,no2gnode[s][b][nd]->cov);

		// 					for(int nd=0;nd<graphno[s][b];nd++) {
		// 						// fprintf(stderr,"Node %d with parents:",i);
		// 						for(int c=0;c<no2gnode[s][b][nd]->child.Count();c++) {
		// 							fprintf(uinigraph_out,"%d->",nd);			
		// 							fprintf(uinigraph_out,"%d;",no2gnode[s][b][nd]->child[c]);
		// 						}
		// 					}
		// 				}
		// 				fprintf(uinigraph_out,"}\n");
		// 				g_idx += 1;
		// 				fprintf(stderr,"g_idx: %d\n", g_idx);
		// 			}
		// 		}
		// 	}
		// }
		/****************
		 **  END KH Adding 
		****************/
        /*****************************
		 ** 3-2. Creating the Unispg for the universal graph and add it into UnispgGp
		*****************************/

        if (multiMode) {
            for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> thoses shouldn't have junctions
                int s=sno/2; // adjusted strand due to ignoring neutral strand				
				unispg_gp->AddGraph(fidx, s, no2gnode[s], bundle[sno].Count());
            }
        }

		/*****************************
		 ** 4. I can clean up some data here:
		*****************************/
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
    				process_transfrags_unispg(s,graphno[s][b],edgeno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],gpos[s][b],guidetrf,pred,trflong);
    				//get_trf_long(graphno[s][b],edgeno[s][b], gpos[s][b],no2gnode[s][b],transfrag[s][b],geneno,s,pred,trflong);


    				// /*
    				{ //DEBUG ONLY
    					//printTime(stderr);
    					fprintf(stderr,">>> There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
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
    					fprintf(stderr,">>> There are %d transfrags[%d][%d]:\n",transfrag[s][b].Count(),s,b);
    					for(int t=0;t<transfrag[s][b].Count();t++) {
    						fprintf(stderr,"%d: ",t);
    						//printBitVec(transfrag[s][b][t]->pattern);
    						fprintf(stderr," %f(%f) long=%d short=%d nodes=%d",transfrag[s][b][t]->abundance,transfrag[s][b][t]->srabund, transfrag[s][b][t]->longread,transfrag[s][b][t]->shortread,transfrag[s][b][t]->nodes.Count());
    						for(int i=0;i<transfrag[s][b][t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[s][b][t]->nodes[i]);
    						if(!transfrag[s][b][t]->abundance) fprintf(stderr," *");
    						fprintf(stderr,"\n");
    					}

    					fprintf(stderr,"\n\n");
    				}
    				// */

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




int infer_transcripts_unispg(BundleData* bundle, int fidx) {
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
	// else if (multiMode && (bundle->keepguides.Count() || !eonly)) {
	// 	fprintf(stderr, "This is the multiMiode of the graph.\n");
	// 	count_good_junctions(bundle);
	// 	geneno = build_graphs_multi(bundle, unispg_gp);
	// }
	else if(bundle->keepguides.Count() || !eonly) {
		//fprintf(stderr,"Process %d reads from %lu.\n",bundle->readlist.Count(),bundle->numreads);

		count_good_junctions(bundle);

		// geneno = build_graphs(bundle);
		geneno = build_graphs_unispg(bundle, fidx);
		// geneno = build_graphs_unispg(bundle, unispg_gp);
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