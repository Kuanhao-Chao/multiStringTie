#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#include "rlink_C.h"


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
		int &edgeno,int &lastgpos,GArray<GEdge>& guideedge, int refend){

	if (multiMode) {
        // int uni_refstart = unispg_gp -> get_refstart();
		// int uni_refend = unispg_gp -> get_refend();
        // fprintf(stderr, "* uni_refstart: %d\n", uni_refstart);
        // fprintf(stderr, "* uni_refend: %d\n", uni_refend);
	}else {
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
		uint endbundle=bundlenode->end; // initialize end with bundlenode's end for now
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
			if(mixedMode) {
				CTransfrag *longtr=new CTransfrag(nodes,trpat,futuretr[i+2]);
				longtr->longread=true;
				transfrag[s][g].Add(longtr);
			}
			else {
				if(longreads)// || mixedMode)
					transfrag[s][g].Last()->longread=true;
			}
		}
	}


	/*****************************
	 ** Step 7: 'traverse_dfs_unispg' function
	 ** 	finished reading bundle -> now create the parents' and children's patterns
	 *****************************/
	GVec<bool> visit;
	visit.Resize(graphno);
	GBitVec parents(graphno+edgeno);

	fprintf(stderr,"traverse graph[%d][%d] now with %d nodes, %d edges and lastgpos=%d....\n",s,g,graphno,edgeno,lastgpos);//edgeno=0;
	traverse_dfs_unispg(s,g,source,sink,parents,graphno,visit,no2gnode,transfrag,edgeno,gpos,lastgpos);
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












/****************************************
 **  1. Add reference pattern to the transfrags as backbone
 **  2. Eliminate transfrags below threshold
 **  3. Add source / sink links (only if they need to be added to explain the traversals)
 **  4. Add edges to disconnected parent-child nodes.
 **  5. For all nodes, check if there is a connection to child.
 **  6. Sort transfrags
 **  7. Create compatibilities
 **  8. 
	 **  set source-to-child transfrag abundances: optional in order not to keep these abundances too low:
	 **  update the abundances of the transfrags coming in from source and going to a node that doesn't have other parents than source
	 **  * this part was removed to improve performance
 ****************************************/
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

	/****************************************
	 ** Step 3:
	 ** add source/sink links but only if they need to be added to explain the traversals in the graph
	 ****************************************/
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

		// /*
		{ // DEBUG ONLY
			//printTime(stderr);
			fprintf(stderr,"\nThere are %d longsrtransfrags after clean up:\n",srfrag.Count());
			for(int i=0;i<srfrag.Count();i++) {
				fprintf(stderr,"srfrag[%d] longstart=%d longend=%d abund=%f t=%d:",i,srfrag[i]->longstart,srfrag[i]->longend,srfrag[i]->abundance,int(srfrag[i]->usepath));
				for(int j=0;j<srfrag[i]->nodes.Count();j++) fprintf(stderr," %d",srfrag[i]->nodes[j]);
				fprintf(stderr,"\n");
			}
		}
		// */

		for(int t1=0;t1<srfrag.Count();t1++) {
			// /*
			fprintf(stderr,"Consider t=%d with abund=%f and nodes:",t1,srfrag[t1]->abundance);
			for(int j=0;j<srfrag[t1]->nodes.Count();j++) {
				if(j) {
					int *pos=gpos[edge(srfrag[t1]->nodes[j-1],srfrag[t1]->nodes[j],gno)];
					if(pos && srfrag[t1]->pattern[*pos]) {
						fprintf(stderr,"-");
					}
					else fprintf(stderr," ");
				}
				fprintf(stderr,"%d",srfrag[t1]->nodes[j]);
			} fprintf(stderr,"\n");
			// */
			if(!srfrag[t1]->nodes[0]) {
				hassource[srfrag[t1]->nodes[1]]=t1;
				fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",srfrag[t1]->nodes[1],t1,srfrag[t1]->abundance);
			}
			else if(srfrag[t1]->nodes.Last()==gno-1) {
				hassink[srfrag[t1]->nodes[0]]=t1;
				fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",srfrag[t1]->nodes[0],t1,srfrag[t1]->abundance);
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
					fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
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
										fprintf(stderr,"trf %d includes %d\n",t[0],t[1]);
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
									fprintf(stderr,"trf %d is intronic including %d\n",t[1],t[0]);
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
								fprintf(stderr,"trf %d %d equivalent start/ends\n",t[1],t[0]);
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
						fprintf(stderr,"keep srfrag %d\n",t1);
					}
					else { // incomplete transcript, possibly wrong
						srfrag[t1]->weak=1;
						fprintf(stderr,"Incomplete transcript %d\n",t1);
					}
				}
			}
		}


		//GBitVec guidesource(gno);
		//GBitVec guidesink(gno);
		//for(int i=0;i<keeptrf.Count();i++) {
		for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier
			fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);
			int n1=srfrag[keeptrf[i].t]->nodes[0];
			int n2=srfrag[keeptrf[i].t]->nodes.Last();
			//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);

			if(!srfrag[keeptrf[i].t]->guide && ((!no2gnode[n1]->hardstart && (hassource[n1]<0 || !keepsource[n1])) ||
					(!no2gnode[n2]->hardend && (hassink[n2]<0 || !keepsink[n2])))) {
			//if((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))) {
				srfrag[keeptrf[i].t]->usepath=-2-ntrf; // order of transfrag
				// fprintf(stderr,"keeptrf[%d].t=%d srfrag.usepath=%d\n",i,keeptrf[i].t,t);
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
					fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
				}
				if(n2==srfrag[keeptrf[i].group[j]]->nodes.Last()) {
					addsink[n2]+=srfrag[keeptrf[i].group[j]]->abundance;
					fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
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

		// /*
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
		// */


		// add source/sink connections
		for(int i=1;i<gno-1;i++) {
			fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f hardstart=%d hardend=%d\n",i,hassource[i],addsource[i],hassink[i],addsink[i],no2gnode[i]->hardstart,no2gnode[i]->hardend);
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
					fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
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
	 ** Step 4:
	 ** add edges between disconnected parent-child nodes
	 ****************************************/
	for(int t=0;t<transfrag.Count();t++) allpat=allpat | transfrag[t]->pattern;

	/****************************************
	 ** Step 5:
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
	 ** Step 6:
	 ** sort transfrag with smallest being the one that has the most nodes, and ties are decided by the abundance (largest abundance first); last transfrags all have 1 node
	 ****************************************/
	if(trsort)
		transfrag.Sort(trCmp);

	// /*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags that remained\n",transfrag.Count());
	}
	// */

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
	 ** Step 7:
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
	 ** Step 8:
	 **  set source-to-child transfrag abundances: optional in order not to keep these abundances too low:
	 **  update the abundances of the transfrags coming in from source and going to a node that doesn't have other parents than source
	 **  * this part was removed to improve performance
	 ****************************************/
	CGraphnode *source=no2gnode[0];
	for(int i=0;i<source->child.Count();i++) {
		float abundance=0;
		int t0=-1;
		if(no2gnode[source->child[i]]->parent.Count()==1 && !no2gnode[source->child[i]]->parent[0]) { // source is the only parent of node
			for(int j=0;j<no2gnode[source->child[i]]->trf.Count();j++) { // iterating all transfrags
				int t=no2gnode[source->child[i]]->trf[j];
				if(transfrag[t]->nodes.Last()==source->child[i]) t0=t; // the current node is the last node of the transfrag
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


/****************************************
 ** 1. process in and out coverages for each node
 ** 2. process guides first
 ** 3. Run maxflow algorithm
 ****************************************/
int find_transcripts_unispg(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand,
		GVec<CGuide>& guidetrf,GPVec<GffObj>& guides,GVec<int>& guidepred,BundleData* bdata,GVec<int>& trflong) {

	GList<CPrediction>& pred = bdata->pred;
	/*
	if(trflong.Count()) get_trf_long_unispg(gno,edgeno, gpos,no2gnode,transfrag,geneno,strand,pred,trflong,bdata);
	if(longreads) return(geneno);*/

	if(longreads) {
		if(trflong.Count()) get_trf_long_unispg(gno,edgeno, gpos,no2gnode,transfrag,geneno,strand,pred,trflong,bdata);
		return(geneno);
	}

	/****************************************
	 ** Step 1:
	 ** process in and out coverages for each node
	 ****************************************/
	int maxi=0; // node with maximum coverage
	GVec<float> nodecov; // node coverages

	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0);
		if(i) { // for all nodes but the source
		    if(i<gno-1 && inode->len()) nodecov[i]=inode->cov/inode->len(); // sink also has 0 coverage
		    if(nodecov[i]>nodecov[maxi]) maxi=i;
		    int nn=inode->trf.Count();
		    float abundin=0;
		    float abundout=0;
		    float abundthrough=0;
		    for(int j=0;j<nn;j++){
		    	int t=inode->trf[j];
		    	if(transfrag[t]->nodes.Last()==i) { // transfrag ends at this node (in transfrag)
		    		abundin+=transfrag[t]->abundance;
		    	}
		    	else if(transfrag[t]->nodes[0]==i) { // transfrag starts at this node (out transfrag)
		    		abundout+=transfrag[t]->abundance;
		    	}
		    	else if(transfrag[t]->pattern[i]) { // through transfrag (here I checked that the transfrag clearly goes through the node)
		    		abundthrough+=transfrag[t]->abundance;
		    	}
		    }
		    if(abundin) inode->rate=abundout/abundin;
		    if(abundout) inode->capacity=abundout+abundthrough; // node capacity tells me how much of that node coverage I can use given how many transfrags leave the node
		    else inode->capacity=abundin+abundthrough;
		} // end if i

		/*
		{ // DEBUG ONLY
			printTime(stderr);
			fprintf(stderr,"Node %d: cov=%f capacity=%f rate=%f ",i,inode->cov/(inode->end-inode->start+1),inode->capacity,inode->rate);
			fprintf(stderr,"trf=");
			for(int t=0;t<inode->trf.Count();t++) fprintf(stderr," %d(%f)",inode->trf[t],transfrag[inode->trf[t]]->abundance);
			fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
		}
		*/
	} // end for i

	GBitVec istranscript(transfrag.Count());
	GBitVec pathpat(gno+edgeno);

/*
#ifdef GMEMTRACE
	double vm,rsm;
	get_mem_usage(vm, rsm);
	GMessage("\t\tM(after istranscript and pathpat init):find_transcripts memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

	/****************************************
	 ** Step 2:
	 ** process guides first
	 ****************************************/
	//fprintf(stderr,"guidetrf.count=%d\n",guidetrf.Count());
	//if(guidetrf.Count()) maxi=guides_flow(gno,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat);
	bool first=true;

	//fprintf(stderr,"guide count=%d\n",guidetrf.Count());

	/****************************************
	 ** Step 3:
	 ** Run maxflow algorithm
	 ****************************************/
	if (eonly) {
		guides_pushmaxflow_unispg(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
	} else { //if(!eonly) {
		if(mixedMode && trflong.Count()) get_trf_long_mix_unispg(gno,edgeno, gpos,no2gnode,transfrag,geneno,strand,pred,trflong,nodecov,istranscript,pathpat,bdata,first);

		if(!mixedMode && guidetrf.Count()) maxi=guides_pushmaxflow_unispg(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
		// /*
		{ // DEBUG ONLY
			if(mixedMode) {
				fprintf(stderr,"After get_trf_long_unispg:\n");
				for(int i=0;i<gno;i++) {
					CGraphnode *inode=no2gnode[i];
					printTime(stderr);
					fprintf(stderr,"Node %d: cov=%f capacity=%f rate=%f ",i,inode->cov/(inode->end-inode->start+1),inode->capacity,inode->rate);
					fprintf(stderr,"trf=");
					for(int t=0;t<inode->trf.Count();t++) fprintf(stderr," %d(%f)",inode->trf[t],transfrag[inode->trf[t]]->abundance);
					fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				}
				fprintf(stderr,"There are %d transfrags:\n",transfrag.Count());
				for(int t=0;t<transfrag.Count();t++) {
					fprintf(stderr,"%d: ",t);
					//printBitVec(transfrag[s][b][t]->pattern);
					fprintf(stderr," %f(%f,%d,%d) long=%d nodes=%d",transfrag[t]->abundance,transfrag[t]->srabund, transfrag[t]->longstart,transfrag[t]->longend,transfrag[t]->longread,transfrag[t]->nodes.Count());
					for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
					if(!transfrag[t]->abundance) fprintf(stderr," *");
					fprintf(stderr,"\n");
				}
			}
		}
		// */
		if(nodecov[maxi]>=1) { // sensitive mode only; otherwise >=readthr
			// 1:
			// parse_trf_weight_max_flow(gno,no2gnode,transfrag,geneno,strand,pred,nodecov,pathpat);
			// 2:
			GBitVec usednode(gno+edgeno);
			parse_trf_unispg(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,usednode,0,pathpat);
		}
	}
	return(geneno);
}
