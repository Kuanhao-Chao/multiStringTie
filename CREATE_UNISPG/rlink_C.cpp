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
			
			transfrag[s][g].Add(tr);
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
 
			transfrag[s][g].Add(tr);
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
	/****************
	 **  KH Adding 
	****************/
	// fprintf(stdout, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	// fprintf(stdout, "&&&&&& Start 'create_graph_unispg'\n");
	// fprintf(stdout, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	/****************
	 **  END KH Adding 
	****************/

	GVec<float>* bpcov = bdata ? bdata->bpcov : NULL; // I might want to use a different type of data for bpcov to save memory in the case of very long bundles

	CGraphnode* source=new CGraphnode(0,0,0);
	no2gnode[s][g].Add(source);
	CGraphnode* sink=new CGraphnode();
	int njunctions=junction.Count();

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

	// This is for reference-guide assembly too.
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

	int graphno=1; // number of nodes in graph (source already included.)
    GIntHash< GVec<int>* > ends; // keeps ids of all nodes ending at a certain position; OR ALL NODES THAT ARE LINKED BY JUNCTIONS TO A CERTAIN POSITION
	GVec<float> futuretr; //future transfrags

	/*****************************
	 ** Step 2: I have a bunch of junctions at the start for which I need to create ends
	 **  This is for mergeMode. Skip here.
	 *****************************/

	/*****************************
	 ** Step 3: 'create_graphnode_unispg' function
	 ** 	Process nodes in the bundle.
	 *****************************/
	int f=0; // feature index
	uint bundle_start=bundlenode->start;
	uint bundle_end=bnode[bundle->lastnodeid]->end;
	// GHashMap<int, int> global2local_nodehash(false); //hash of pointers
	// global2local_nodehash.Add(0, 0);
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
			if((nje<njunctions && ejunction[nje]->end - currentstart < junctionsupport) &&
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
					// fprintf(stderr,"1 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
				}
			}
			else { // I haven't seen nodes before that finish here (maybe due to error correction?) => link to source
				source->child.Add(graphnode->nodeid);  // this node is the child of source
				graphnode->parent.Add(source->nodeid); // this node has source as parent
				// COUNT EDGE HERE
				edgeno++;
				// fprintf(stderr,"2 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
			}
		}
		else { // this node comes from source directly
			source->child.Add(graphnode->nodeid);  // this node is the child of source
			graphnode->parent.Add(source->nodeid); // this node has source as parent
			// COUNT EDGE HERE
			edgeno++;
			// fprintf(stderr,"3 Edge 0-%d, edgeno=%d\n",graphnode->nodeid,edgeno);
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
			// if(nje<njunctions) fprintf(stderr,"Found junction:%d-%d(%d)\n",ejunction[nje]->start,ejunction[nje]->end,ejunction[nje]->strand);


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
				if(trim) graphnode=trimnode_all(s,g,refstart,junction[njs]->start,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


				dropcov=true;
				// if no trimming required just set the end of the node
				graphnode->end=junction[njs]->start; // set the end of current graphnode to here; introduce smaller nodes if trimming is activated
				uint pos=junction[njs]->start;
				while(njs<njunctions && junction[njs]->start==pos ) { // remember ends here
					if((junction[njs]->strand+1) == 2*s) {
						//seenjunc++;
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
					njs++;
				}

				if(pos<endbundle) { // there is still place for another node in this bundle (I might put a limit of length here for the graphnode -> because otherwise one can assume this is just a pre-mRNA fragment)
					// see if I should just skip node
					if(endbundle-pos<junctionsupport) {
						while(njs<njunctions && junction[njs]->strand+1 != 2*s) njs++;
						if((njs>=njunctions || junction[njs]->start > endbundle) && (nje>=njunctions || ejunction[nje]->end > endbundle)) { // there are no more junctions starting within this bundle
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
						// fprintf(stderr,"4 Edge %d-%d, edgeno=%d nextnode: %u-%u pos=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno,nextnode->start,nextnode->end,pos);
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

								if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
								else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

							}
						}
					}

					if(trim) graphnode=trimnode_all(s,g,refstart,pos-1,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);// do something to find intermediate nodes; alternatively, I could only do this for end nodes


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
					// fprintf(stderr,"5 Edge %d-%d, edgeno=%d\n",graphnode->nodeid,nextnode->nodeid,edgeno);

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
					// fprintf(stderr,"6 Edge %d-%d, edgeno=%d\n",node->nodeid,graphnode->nodeid,edgeno);
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

						if(sourceguide)	graphnode=source2guide(s,g,refstart,start,end,graphnode,source,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);
						else graphnode=guide2sink(s,g,refstart,start,end,graphnode,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno);

					}
				}
			}

			if(trim) graphnode=trimnode_all(s,g,refstart,endbundle,graphnode,source,sink,bpcov,futuretr,graphno,bundlenode,bundle2graph,no2gnode,edgeno); // do something to find intermediate nodes; alternatively, I could only do this for end nodes


			graphnode->end=endbundle;
			// COUNT EDGE HERE (this is an edge to sink)
			edgeno++;
			// fprintf(stderr,"7 Edge to sink from %d, edgeno=%d\n",graphnode->nodeid,edgeno);
		}
	    bundlenode=bundlenode->nextnode; // advance to next bundle
	} // end while(bundlenode!=NULL)
    

	/*****************************
	 ** Step 4: 'get_cov_sign' function
	 ** 	add source/sink links for very high coverage drops
	 *****************************/
	for(int i=1;i<graphno;i++) {
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
	if(graphno>allowed_nodes) { // TODO: define allowed_nodes as a default in stringtie.cpp that varies with the memory
		graphno=prune_graph_nodes(graphno,s,g,bundle2graph,bnode.Count(),no2gnode,junction,edgeno,futuretr,sink);
	}

	sink->nodeid=graphno;
	no2gnode[s][g].Add(sink);
	graphno++;

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
			transfrag[s][g].Add(tr);
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
	/*
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
	*/
	/****************************************
	 ** Step 1:
	 ** add all guide patterns to the set of transfrags so that I can have a "backbone" for each guide
	 ** I need this because there might be an incompatible transfrag connecting the nodes in the guide 
	 ****************************************/
	//fprintf(stderr,"There are %d guides\n",guidetrf.Count());
	for(int i=0;i<guidetrf.Count();i++) if(guidetrf[i].trf->guide){

		CTransfrag *t=NULL;
		bool add=true;
		t=new CTransfrag(guidetrf[i].trf->nodes,guidetrf[i].trf->pattern,trthr*ERROR_PERC);

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


		t->guide=1+guidetrf[i].g;
		//fprintf(stderr,"t->guide set to=%d\n",1+guidetrf[i].g);
		t->longstart=no2gnode[t->nodes[0]]->start;
		t->longend=no2gnode[t->nodes.Last()]->end;
		no2gnode[t->nodes[0]]->hardstart=1;  // I can always trust a guide's start
		no2gnode[t->nodes.Last()]->hardend=1; // I can always trust a guide's end
		if(add) transfrag.Add(t);
	}
	else {
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

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d] abund=%f long:%d: short:%d usepath:%d",i,transfrag[i]->abundance,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	*/

	GBitVec allpat(gno+edgeno);

	bool trsort=true;

	/****************************************
	 ** Step 3:
	 ** add source/sink links but only if they need to be added to explain the traversals in the graph
	 ****************************************/

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
	 ** Step 7:
	 ** create compatibilities
	 ****************************************/
	for(int t1=0;t1<transfrag.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations
		// update nodes
		int n1=transfrag[t1]->nodes.Count();

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
		/*
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}
		*/
	} // end for(int t1=0;t1<transfrag.Count();t1++)



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


	/****************************************
	 ** Step 3:
	 ** Run maxflow algorithm
	 ****************************************/
	if (eonly) {
		guides_pushmaxflow_unispg(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
	} else { //if(!eonly) {
		
		if(guidetrf.Count()) maxi=guides_pushmaxflow_unispg(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
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
