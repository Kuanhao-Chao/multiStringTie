// #include "rlink.h"
#include "rlink_C_help.h"

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

void parse_trf_unispg(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& usednode,float maxcov,GBitVec& prevpath) {

	 GVec<int> path;
	 path.Add(maxi);
	 GBitVec pathpat(gno+edgeno);
	 pathpat[maxi]=1;
	 istranscript.reset();

	 float flux=0;
	 //float fragno=0;
	 GVec<float> nodeflux;

	 /*
	 { // DEBUG ONLY
	 	 fprintf(stderr,"\n\n***Start parse_trf_unispg with maxi=%d and cov=%f\n",maxi,nodecov[maxi]);
		 //fprintf(stderr,"Transcripts before path:");
		 //for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		 //fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf_unispg memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 }
	 */


	if(back_to_source_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
		 	 if(includesource) path.cAdd(0);
	 		 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			if(fwd_to_sink_fast(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
	 			 bool full=true;

				 flux=push_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,gpos,full);

				 /*
	 			 { // DEBUG ONLY
	 				 //printTime(stderr);
	 				 fprintf(stderr,"flux=%g Path:",flux);
	 				 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
	 				 fprintf(stderr,"***\n");
	 			 }
	 			 */
	 		}
			/*else {
	 			//pathpat.reset();
	 			//pathpat[maxi]=1;
	 		}*/
	 }
	 /*else {
		 //pathpat.reset();
		 //pathpat[maxi]=1;
	 }*/

	 //bool cont=true;

	 if(flux>epsilon) {
		 bool included=true;
		 float cov=store_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,strand,gno,gpos,included,prevpath);

		 /*
		 { // DEBUG ONLY
			 //fprintf(stderr,"Prevpath=");
			 //printBitVec(prevpath);
			 //fprintf(stderr,"\n");
		 	 fprintf(stderr,"cov=%f maxcov=%f\n",cov,maxcov);
		 }
		 */

		 float frac=isofrac;
		 if(mixedMode && frac<ERROR_PERC) frac=ERROR_PERC;

		 if(included || cov<frac*maxcov) {
			 /*
			 if(sensitivitylevel) usednode[maxi]=1;
			 else usednode = usednode | prevpath;
			 */
			 //fprintf(stderr,"included\n");
			 usednode[maxi]=1;
			 //fprintf(stderr,"used maxi=%d\n",maxi);
			 maxi=0;
			 maxcov=0;
			 //cont=false;

		 }
		 else if(cov>maxcov) maxcov=cov;

	 }
	 else {
		 /*
		 if(sensitivitylevel) usednode[maxi]=1; // start at different locations in graph
		 else {
			 usednode = usednode | prevpath;
			 usednode = usednode | pathpat;
		 }
		 */
		 usednode[maxi]=1;

		 //fprintf(stderr,"set usednode of %d\n",maxi);

		 maxi=0;
		 maxcov=0;
		 //cont=false;
	 }

	 //fprintf(stderr,"maxcov=%g gno=%d\n",maxcov,gno);

	 // Node coverages:
	 for(int i=1;i<gno;i++)
		 if(!usednode[i] && nodecov[i]>nodecov[maxi]) maxi=i;

	 //fprintf(stderr," maxi=%d nodecov=%f\n",maxi,nodecov[maxi]);

	 //if(nodecov[maxi]>=readthr && (!specific || cont)) { // if I still have nodes that are above coverage threshold
	 if(nodecov[maxi]>=1) { // if I still have nodes to use

		 /*
		 { // DEBUG ONLY
			 //printTime(stderr);
			 fprintf(stderr,"\nAfter update:\n");
			 for(int i=0;i<gno;i++) {
				 fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
				 fprintf(stderr,"trf=");
				 for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
				 fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
			 }
		 }
		 */

		 path.Clear();
		 nodeflux.Clear();
		 parse_trf_unispg(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,usednode,maxcov,prevpath);
	 }
}

int guides_pushmaxflow_unispg(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata) {

	int maxi=1;
	int ng=guidetrf.Count();

	if(ng==1) { // if only one guide I do not need to do the 2 pass
		GVec<float> nodeflux;
		//float fragno=0;
		bool full=true;
		float flux= push_max_flow(gno,guidetrf[0].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[0].trf->pattern,gpos,full);
		istranscript.reset();

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"guide=%s flux[0]=%g\n",guides[guidetrf[0].g]->getID(),flux);
		}
		*/

		if(flux>epsilon) {
			bool include=true;
			if(guidepred[guidetrf[0].g]==-1) {

				store_transcript(pred,guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[0].g]);
				//if(eonly) { // this is not correct because it might have been assigned before
					guidepred[guidetrf[0].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
					//fprintf(stderr,"guidepred[%d]=%d\n",guidetrf[0].g,guidepred[guidetrf[0].g]);
				//}
			}
			else {
				update_guide_pred(pred,guidepred[guidetrf[0].g],guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
			}

		}

		//if(nodecov[maxi]<readthr) break; // no need to find other paths since they aren't any above allowed read threshold
		//if(nodecov[maxi]<1) break; // I shouldn't be restricting this at all?


		/*
		{ // DEBUG ONLY
		  fprintf(stderr,"\nAfter update:\n");
		  for(int i=0;i<gno;i++) {
			  fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
			  fprintf(stderr,"trf=");
			  for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
			  fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
		  }
		}
		*/

		//return(maxi);
	}
	else if(ng) {


		// this is an abundance for guides based on maximum flow for each guide (how much can they each carry
		for(int g=0;g<guidetrf.Count();g++) {
			guidetrf[g].trf->abundance=push_guide_maxflow(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,guidetrf[g].trf->pattern);
			istranscript.reset();
		}
		guidetrf.Sort(guidedabundCmp);

		/*
		{ // DEBUG ONLY
			for(int g=0;g<guidetrf.Count();g++) {
				fprintf(stderr,"Abundance of guide[%d]=%f with nodes:",g,guidetrf[g].trf->abundance);
				for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) fprintf(stderr," %d",guidetrf[g].trf->nodes[i]);
				fprintf(stderr,"\n");
			}
		}
		*/

		GVec<float> nodeflux;
		for(int g=ng-1;g>=0;g--) if(guidetrf[g].trf->abundance){ // calculate maximum push flow for each guide starting from the less covered one

			float flux=guidepushflow(g,guidetrf,gno,istranscript,transfrag,no2gnode,nodeflux);

			istranscript.reset();

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"guide=%s flux[%d]=%f\n",guides[g]->getID(),g,flux);
			}
			*/

			bool include=true;
			if(flux>epsilon) {
				if(guidepred[guidetrf[g].g]==-1) {
					store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[g].g]);
					//if(eonly) {
						guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
						//fprintf(stderr,"2 guidepred[%d]=%d\n",guidetrf[g].g,guidepred[guidetrf[g].g]);
						//}

						/*
						{ // DEBUG ONLY
				  	  	  fprintf(stderr,"\nAfter update:\n");
				  	  	  for(int i=0;i<gno;i++) {
					  	  	  fprintf(stderr,"Node %d: %f ",i,nodecov[i]);
					  	  	  fprintf(stderr,"trf=");
					  	  	  for(int t=0;t<no2gnode[i]->trf.Count();t++) fprintf(stderr," %d(%f)",no2gnode[i]->trf[t],transfrag[no2gnode[i]->trf[t]]->abundance);
					  	  	  fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
				  	  	  }
						}
						 */
				}
				else {
					update_guide_pred(pred,guidepred[guidetrf[g].g],guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
				}

				nodeflux.Clear();

			}
			else { // it's possible that this is a single exon gene included in a much larger interval and this is why it didn't get predicted
				if(guides[guidetrf[g].g]->exons.Count()==1) { // single exon gene included in another prediction
					// check if it overlaps other single exon genes

					bool overlap=false;
					for(int r=ng-1;r>g;r--) if(guidepred[guidetrf[r].g]>-1){
						if((guidetrf[g].trf->pattern & guidetrf[r].trf->pattern)==guidetrf[g].trf->pattern) {
							if(guides[guidetrf[r].g]->exons.Count()>1) {
								overlap=true;
								break;
							}
							else { // overlaps single exon gene -> check if it's a true overlap
								if(guides[guidetrf[r].g]->exons[0]->overlap(guides[guidetrf[g].g]->exons[0])) {
									overlap=true;
									break;
								}
							}
						}
					}
					if(!overlap) { // no overlap detected
						for(int z=0;z<guidetrf[g].trf->nodes.Count();z++) nodeflux.cAdd(1.0);
						if(guidepred[guidetrf[g].g]==-1) {
							store_transcript(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[g].g]);
							//if(eonly) {
								guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript
								//fprintf(stderr,"2 guidepred[%d]=%d\n",guidetrf[g].g,guidepred[guidetrf[g].g]);
							//}
						}
						else {
							update_guide_pred(pred,guidepred[guidetrf[g].g],guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
						}
						nodeflux.Clear();
					}
				}
			}
		}

	}

	// Node coverages:
	for(int i=1;i<gno-1;i++)
		if(nodecov[i]>nodecov[maxi]) maxi=i;


	if(eonly && nodecov[maxi]>epsilon) { // this is the end for eonly so I should make sure I use all reads

		for(int g=0;g<guidetrf.Count();g++) delete guidetrf[g].trf;
		guidetrf.Clear();

		char strand='-';
		if(s) strand='+';

		GVec<CNodeGuide> nodeinfo(gno);
		nodeinfo.setCount(gno);

		int ng=0;

		GVec<int> olen;
		// find guides' partial patterns
		for(int g=0;g<guides.Count();g++) {
			//fprintf(stderr,"Consider guide[%d out of %d] %s\n",g,guides.Count(),guides[g]->getID());
			if((guides[g]->strand==strand) && (guides[g]->overlap(no2gnode[1]->start,no2gnode[gno-2]->end))) { // if guide on the same strand and overlaps at all the graph
				//fprintf(stderr,"...there is overlap\n");
				olen.Clear();
				int olensum;
				CTransfrag *trguide=find_guide_partial_pat(guides[g],no2gnode,gno,edgeno,gpos,olen,olensum);
				if(trguide) { // the guide can be found among the graph nodes
					//fprintf(stderr,"...partial pattern found!\n");
					CGuide newguide(trguide,g);
					guidetrf.Add(newguide);
					for(int i=0;i<trguide->nodes.Count();i++) {
						CPartGuide pg(ng,olen[i],olensum,guides[g]->covlen);
						nodeinfo[trguide->nodes[i]].guide.Add(pg);
						int idx=nodeinfo[trguide->nodes[i]].guide.Count()-1; // position of guide in nodeinfo
						bool terminal=true;
						if(i) { // not first node in guide covered
							int *pos=gpos[edge(trguide->nodes[i-1],trguide->nodes[i],gno)];
							if(pos && trguide->pattern[*pos]) terminal=false; // if there is edge from previous node and it's present in guide
						}
						nodeinfo[trguide->nodes[i]].guide[idx].terminal_in=terminal;
						terminal=true;
						if(i<trguide->nodes.Count()-1) { // not last node in guide covered
							int *pos=gpos[edge(trguide->nodes[i],trguide->nodes[i+1],gno)];
							if(pos && trguide->pattern[*pos]) terminal=false; // if there is edge from previous node and it's present in guide
						}
						nodeinfo[trguide->nodes[i]].guide[idx].terminal_out=terminal;
					}
					ng++;
				}
			}
		}

		// future: compute maxflow for each guide starting at the most covered node maybe and going all the way to where I can do it.
		// or maybe compute maxflow for each guide from all terminal nodes and then assign coverages based on winners or the one that gets the most coverage
		// hint: max_flow doesn't need source/sink in the path so I can go with that
		// or: do some type EM algorithm like I am trying to do below

		GVec<int> path(1);
		path.setCount(1);
		GVec<float> nodeflux(1);
		nodeflux.setCount(1);

		GVec<int> coverednode; // array of all nodes that still have coverage _and_ guides

		// determine gcounts and the strict assignments to guides
		for(int n=1;n<gno-1;n++) if(nodecov[n]>epsilon && nodeinfo[n].guide.Count()) { // for all nodes for which there is still coverage and have guides
			ng=nodeinfo[n].guide.Count(); // number of guides in node
			if(ng==1) { // there is only one guide -> it gets all the coverage
				int g=guidetrf[nodeinfo[n].guide[0].idx].g;
				path[0]=n;
				nodeflux[0]=1;
				//fprintf(stderr,"1 g=%d\n",g);
				if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
				else // new prediction for this guide
					guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
			}
			else { // there are more than one guide overlapping the node -> count the strict/loose counts; might want to rethink doing this for single nodes in the graph

				float incount=0; // should be distrtibuted between incoming guides that continue and terminal ones -> not doing this now
				float outcount=0;
				//float throughcount=0;
				float terminal_incount=0;
				float terminal_outcount=0;
				float guide_incount=0;
				float guide_outcount=0;

				GVec<float> in_strictcount(ng,float(0)); // these are abundances that are unique to a single guide -> this will get assigned to guides
				GVec<float> in_loosecount(ng,float(0)); // these are guide abundances that are shared with other guides
				GVec<float> in_looseterminalcount(ng,float(0)); // these are guide abundances that are shared with other guides
				GVec<float> out_strictcount(ng,float(0));
				GVec<float> out_loosecount(ng,float(0));
				GVec<float> out_looseterminalcount(ng,float(0));

				GVec<CTrGuidePat> out_trcount; // I can use the node trcount for the in_counts but I need for the outcount because it needs to be adjusted by the rate

				int ncovguides=0;
				bool allstrictzero=true;
				bool allloosezero=true;

				CGraphnode *inode=no2gnode[n];
				int nn=inode->trf.Count();
				for(int j=0;j<nn;j++){
					int t=inode->trf[j];
					if(transfrag[t]->abundance>epsilon) { // only if transfrag still has abundance it's worth considering
						GVec<int> compguide; // keeps all guides that are compatible with this transfrag
						GBitVec guidepat(guidetrf.Count()); // pattern of guides that are present in transfrag
						bool terminal_in=true;  // true only if all guides compatible with transfrag are terminal
						bool terminal_out=true;

						for(int i=0;i<ng;i++) {
							int g=nodeinfo[n].guide[i].idx;
							if(((transfrag[t]->pattern) & guidetrf[g].trf->pattern) == transfrag[t]->pattern) { // transfrag is compatible to guide
								compguide.Add(i); // make sure that later you change it to guidetrf indexes
								guidepat[g]=1;
								if(terminal_in && !nodeinfo[n].guide[i].terminal_in) terminal_in=false;
								if(terminal_out && !nodeinfo[n].guide[i].terminal_out) terminal_out=false;
								if(!nodeinfo[n].guide[i].gcount) {
									nodeinfo[n].guide[i].gcount=1;
									ncovguides++;
								}
							}
						}
						if(!compguide.Count()) {
							terminal_in=false;
							terminal_out=false;
						}

						if(transfrag[t]->nodes.Last()==n) { // transfrag ends at this node (in transfrag)
							if(terminal_in)  {
								if(compguide.Count()>1) { // all guides compatible with transfrag are terminal
									terminal_incount+=transfrag[t]->abundance;
									for(int i=0;i<compguide.Count();i++) {
										in_looseterminalcount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
							else incount+=transfrag[t]->abundance; // this biases in favor of continuing guides -> might want to rethink this

							if(compguide.Count()==1) { // only one guide compatible with transfrag
								in_strictcount[compguide[0]]+=transfrag[t]->abundance;
								allstrictzero=false;
							}
							else if(compguide.Count()) { // there are multiple guides in transfrag
								int p=find_cguidepat(guidepat,nodeinfo[n].trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_in);
									nodeinfo[n].trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) nodeinfo[n].trcount[nodeinfo[n].trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else nodeinfo[n].trcount[p].abund+=transfrag[t]->abundance;

								if(!terminal_in) {
									guide_incount+=transfrag[t]->abundance; // terminals can only be guide transfrags
									for(int i=0;i<compguide.Count();i++) {
										in_loosecount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
						}
						else if(transfrag[t]->nodes[0]==n) { // transfrag starts at this node (out transfrag)
							if(terminal_out) {
								if(compguide.Count()>1) {
									terminal_outcount+=transfrag[t]->abundance;
									for(int i=0;i<compguide.Count();i++) {
										out_looseterminalcount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
							else outcount+=transfrag[t]->abundance;
							if(compguide.Count()==1) {
								out_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								allstrictzero=false;
							}
							else if(compguide.Count()){
								int p=find_cguidepat(guidepat,out_trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_out);
									out_trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) out_trcount[out_trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else out_trcount[p].abund+=transfrag[t]->abundance;

								if(!terminal_out) {
									guide_outcount+=transfrag[t]->abundance;

									for(int i=0;i<compguide.Count();i++) {
										out_loosecount[compguide[i]]+=transfrag[t]->abundance;
										allloosezero=false;
									}
								}
							}
						}
						else if(transfrag[t]->pattern[n]) { // through transfrag (here I checked that the transfrag clearly goes through the node)
							//throughcount+=transfrag[t]->abundance;
							incount+=transfrag[t]->abundance; // I do this instead of the above because the throughout could be very high and the in/outcounts shouldn't dominate in that case
							outcount+=transfrag[t]->abundance;

							if(compguide.Count()==1) {
								in_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								out_strictcount[compguide[0]]+=transfrag[t]->abundance; // only one guide compatible with transfrag
								allstrictzero=false;
							}
							else if(compguide.Count()) {
								guide_incount+=transfrag[t]->abundance;
								guide_outcount+=transfrag[t]->abundance;
								int p=find_cguidepat(guidepat,nodeinfo[n].trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_in);
									nodeinfo[n].trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) nodeinfo[n].trcount[nodeinfo[n].trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else {
									nodeinfo[n].trcount[p].abund+=transfrag[t]->abundance;
								}

								p=find_cguidepat(guidepat,out_trcount);
								if(p==-1) { // didn't find the guide pattern
									CTrGuidePat pat(guidepat,transfrag[t]->abundance,terminal_out);
									out_trcount.Add(pat);
									for(int i=0;i<compguide.Count();i++) out_trcount[out_trcount.Count()-1].g.Add(nodeinfo[n].guide[compguide[i]].idx);
								}
								else {
									out_trcount[p].abund+=transfrag[t]->abundance;
								}

								for(int i=0;i<compguide.Count();i++) {
									in_loosecount[compguide[i]]+=transfrag[t]->abundance;
									out_loosecount[compguide[i]]+=transfrag[t]->abundance;
									allloosezero=false;
								}
							}
						}
					}
				} // done computing counts for node

				if(allstrictzero && allloosezero) { // only single guides or no transfrag compatible to guides
												// -> I need to choose some preferences: bias toward longest overlapping guide
					for(int i=0;i<ng;i++) { // ng>1
						nodeinfo[n].guide[i].gcount=1; // give them a uniform probability at node level: coverage probabilities take priority
					}
					nodeinfo[n].sumtrcount=1; // one guide is sufficient to explain all abundances;
					// Note: pattern of guides in trcount should all be with only one 1 for the actual guide but we don't need to use it since using the guide explains everything
				}
				else { // some guides are compatible to the node transfrags

					if(ncovguides==1) { // only one guide is compatible with the node transfrags -> takes all node coverage (this biases strongly against single exon guides, or incomplete guides -> might want to rethink)
						int g=guidetrf[nodeinfo[n].guide[0].idx].g;
						path[0]=n;
						nodeflux[0]=1;
						//fprintf(stderr,"2 g=%d\n",g);
						if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
						else // new prediction for this guide
							guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
					}
					else { // more than one guide got counts

						// first get the rate
						float rate=1;
						nodeinfo[n].sumtrcount=0;
						if(incount && outcount) {
							rate=incount/outcount;
							nodeinfo[n].sumtrcount+=rate*terminal_outcount;
							nodeinfo[n].sumtrcount+=terminal_incount;
						}
						nodeinfo[n].sumtrcount+=guide_incount+rate*guide_outcount+terminal_incount+rate*terminal_outcount;

						float sumstrict=nodeinfo[n].sumtrcount;
						if(!allstrictzero) { // there are strict counts that need to be assigned to
							for(int i=0;i<ng;i++) {
								sumstrict+=in_strictcount[i]+rate*out_strictcount[i];
								if(nodeinfo[n].guide[i].terminal_in) sumstrict+=in_strictcount[i]; // if node is terminal the counts should be adjusted
								if(nodeinfo[n].guide[i].terminal_out) sumstrict+=rate*out_strictcount[i];
							}
						}

						// update guide coverages with unique reads and gcounts for later
						for(int i=0;i<ng;i++) {

							if(!allstrictzero) {
								float strictcount=in_strictcount[i]+rate*out_strictcount[i];
								if(incount && outcount) {
									if(nodeinfo[n].guide[i].terminal_in) strictcount+=in_strictcount[i];
									if(nodeinfo[n].guide[i].terminal_out) strictcount+=rate*out_strictcount[i];
								}
								if(strictcount) { // if there are reads unique to the guide
									int g=guidetrf[nodeinfo[n].guide[i].idx].g;
									path[0]=n;
									nodeflux[0]=strictcount/sumstrict;
									//fprintf(stderr,"3 g=%d\n",g);
									if(guidepred[g]!=-1) update_guide_pred(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,false);
									else // new prediction for this guide
										guidepred[g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],false);
								}
							}

							// now set gcounts
							nodeinfo[n].guide[i].gcount=in_loosecount[i]+rate*out_loosecount[i]+
									in_looseterminalcount[i]+rate*out_looseterminalcount[i];
							if(incount && outcount) {
								if(nodeinfo[n].guide[i].terminal_in) nodeinfo[n].guide[i].gcount+=in_looseterminalcount[i];
								if(nodeinfo[n].guide[i].terminal_out) nodeinfo[n].guide[i].gcount+=rate*out_looseterminalcount[i];
							}
						}

						if(!allstrictzero) { // update node coverage here
							if(allloosezero) nodecov[n]=0;
							else nodecov[n]*=nodeinfo[n].sumtrcount/sumstrict;
						}

						if(!allloosezero && nodecov[n]) { // get the pattern counts
							if(incount && outcount) for(int i=0;i<nodeinfo[n].trcount.Count();i++)
								if(nodeinfo[n].trcount[i].terminal) nodeinfo[n].trcount[i].abund*=2;
							for(int i=0;i<out_trcount.Count();i++) {
								int t=find_cguidepat(out_trcount[i].pat,nodeinfo[n].trcount);
								out_trcount[i].abund*=rate;
								if(out_trcount[i].terminal && incount && outcount) out_trcount[i].abund*=2;
								//else out_trcount[i].terminal=false; // I can not have a trcount with the same guides that is terminal both for in and out
								if(t==-1) { // didn't find the guide pattern
									nodeinfo[n].trcount.Add(out_trcount[i]);
								}
								else nodeinfo[n].trcount[t].abund+=out_trcount[i].abund;
							}
						}

					}
				}
			}
			if(nodecov[n]>epsilon) { // there is still coverage left for the node
				coverednode.Add(n);
			}
		}

		if(coverednode.Count()) { // I still have covered nodes: asign reads to guides; I need to do some sort of EM algorithm here
			// my prior bias in assigning reads should be based on coverages allready found

			int nEM=10; // number of times to repeat the EM algorithm: need to see if this is a good count, or if I should use and epsilon change to stop
			ng=guidetrf.Count();
			GVec<float> initgcov(ng,float(0));

			// first set the initial coverages in covered nodes
			for(int i=0;i<ng;i++) {
				int np=guidepred[guidetrf[i].g];
				if(np!=-1) {
					initgcov[i]=pred[np]->cov*abs(pred[np]->tlen);
				}
			}

			GVec<float> prevgcov(ng,float(0)); // the new gcov that will be computed as part of the EM algorithm

			// now set the cov's in CPartGuide -> for sort purposes and for the EM algorithm
			for(int i=0;i<coverednode.Count();i++) {
				int n=coverednode[i];
				for(int j=0;j<nodeinfo[n].guide.Count();j++) {
					int g=nodeinfo[n].guide[j].idx;
					nodeinfo[n].guide[j].cov=initgcov[g];
					prevgcov[g]=initgcov[g];
				}
			}

			GVec<float> gcov(ng,float(0)); // the new gcov that will be computed as part of the EM algorithm

			int m=0;
			while(m<nEM) { // do the EM

				for(int i=0;i<ng;i++) gcov[i]=initgcov[i];

				for(int i=0;i<coverednode.Count();i++) { // for each node recompute probabilities
					int n=coverednode[i];
					nodeinfo[n].guide.Sort(partguideCmp); // I need to sort at each step because the cov's get updated
					float totalcount=0;
					bool covnotreached=true;
					float sumgcov=0;
					for(int j=0;j<nodeinfo[n].guide.Count();j++) { // for each guide in the node, starting from the most abundant to the least
						nodeinfo[n].guide[j].ncov=0; // I need to set this for each guide -> this is why I not break from the loop below
						int g=nodeinfo[n].guide[j].idx; // index in guidetrf
						sumgcov+=prevgcov[g];
						if(covnotreached) for(int t=0;t<nodeinfo[n].trcount.Count();t++) { // for each trguidepat in the node
							if(nodeinfo[n].trcount[t].pat[g]) { // current guide is in this pattern
								float sum=0; // sum of all coverages of guides in the transcript
								for(int k=0;k<nodeinfo[n].trcount[t].g.Count();k++) {
									int l=nodeinfo[n].trcount[t].g[k];
									sum+=prevgcov[l];
								}
								float abund=nodeinfo[n].sumtrcount-totalcount; // this is the maximum abundance allowed for this guide
								if(sum) { // I have some guides that have coverage
									float newabund=nodeinfo[n].trcount[t].abund*prevgcov[g]/sum;
									if(newabund<abund) abund=newabund;
								}
								else { // no guide has assigned coverages yet -> winner takes all
									if(nodeinfo[n].trcount[t].abund<abund) abund=nodeinfo[n].trcount[t].abund;
								}
								nodeinfo[n].guide[j].ncov+=abund;
								gcov[g]+=abund*nodeinfo[n].guide[j].olen/no2gnode[n]->len();
								totalcount+=abund;
								if(nodeinfo[n].sumtrcount-totalcount<epsilon) {
									covnotreached=false;
									break;
								}
							}
						}
					}
					if(covnotreached) { // there were not enough transfrags to give coverages -> redistribute reads based on coverages
						float abund=nodeinfo[n].sumtrcount-totalcount; // this is how much is left
						if(sumgcov) for(int j=0;j<nodeinfo[n].guide.Count();j++) {
							if(!nodeinfo[n].guide[j].cov) break;
							float newabund=abund*nodeinfo[n].guide[j].cov/sumgcov;
							nodeinfo[n].guide[j].ncov+=newabund;
							int g=nodeinfo[n].guide[j].idx;
							gcov[g]+=newabund*nodeinfo[n].guide[j].olen/no2gnode[n]->len();
						}
						else { // all guides have cov zero -> winner takes all
							nodeinfo[n].guide[0].ncov+=abund;
							int g=nodeinfo[n].guide[0].idx;
							gcov[g]+=abund*nodeinfo[n].guide[0].olen/no2gnode[n]->len();
						}
					}
				} // end node

				bool nochange=true;
				// check coverages
				for(int i=0;i<coverednode.Count();i++) {
					int n=coverednode[i];
					for(int j=0;j<nodeinfo[n].guide.Count();j++) {
						int g=nodeinfo[n].guide[j].idx;
						float diff=nodeinfo[n].guide[j].cov-gcov[g];
						if(diff<0) diff=-diff;
						if(diff) {
							prevgcov[g]=gcov[g];
							nodeinfo[n].guide[j].cov=gcov[g];
							if(diff>1) nochange=false;
						}
					}
				}

				// see if there is a change in probability
				if(nochange) break;

				m++;
			} // end EM

			// new coverages are estimated so now I can assign the reads to the predictions
			for(int g=0;g<guidetrf.Count();g++) if(gcov[g]>initgcov[g]) { // for each guide that has more coverage to add
				path.Clear();
				nodeflux.Clear();
				for(int i=0;i<guidetrf[g].trf->nodes.Count();i++) if(nodecov[guidetrf[g].trf->nodes[i]]) { // for each covered node
					int n=guidetrf[g].trf->nodes[i];
					for(int j=0;j<nodeinfo[n].guide.Count();j++) if(g==nodeinfo[n].guide[j].idx){ // found my guide
						if(nodeinfo[n].guide[j].ncov) { // if there is coverage assigned to node -> i should add it to the path
							path.Add(n);
							float prop=nodeinfo[n].guide[j].ncov/nodeinfo[n].sumtrcount;
							nodeflux.Add(prop);
						}
						break;
					}
				}
				//fprintf(stderr,"4 g=%d\n",guidetrf[g].g);
				if(guidepred[guidetrf[g].g]!=-1) update_guide_pred(pred,guidepred[guidetrf[g].g],path,nodeflux,nodecov,no2gnode,gno,false);
				else // new prediction for this guide
					guidepred[guidetrf[g].g]=store_guide_transcript(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[guidetrf[g].g],false);
			}
		}

		maxi=0;
	}

	return(maxi);
}

void get_trf_long_mix_unispg(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int& geneno,int strand,
		GList<CPrediction>& pred,GVec<int>& trflong,GVec<float>& nodecovall,GBitVec& istranscript,GBitVec& prevpath,BundleData *bdata,bool &first) {

	GPVec<GffObj>& guides = bdata->keepguides;
	GVec<float> nodecov; // the coverage of all transfrags entering a node
	GVec<float> noderate;
	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0); // this are all transfrags that link nodes together
		float rate=1;
		if(i && i<gno-1) {
			for(int j=0;j<inode->trf.Count();j++) { // for all transfrags going through node
				int t=inode->trf[j];
				//if(transfrag[t]->longread && transfrag[t]->nodes[0]<i) { // entering transfrags:
				if(transfrag[t]->longread && transfrag[t]->nodes.Last()>i) { // exiting transfrags: this is more consistent with the nodeflux computation
					nodecov[i]+=transfrag[t]->abundance;
				}
			}

			if(nodecov[i]) rate=nodecov[i];
			if(rate<=0) rate=1; // this shouldn't happen
			//fprintf(stderr,"rate=%f\n",rate);
			rate=inode->cov/rate;
		}
		noderate.Add(rate);
		//fprintf(stderr,"Node[%d]:%d-%d no2gnode->cov=%f nodecov=%f noderate=%f\n",i,inode->start,inode->end,inode->cov,nodecov[i],noderate[i]);
	}

	GVec<int> path;
	GBitVec pathpat(gno+edgeno);
	int minp;
	int maxp;
	int maxi;


	char sign='-';
	if(strand) { sign='+';}
	int npred=pred.Count();

	GVec<CTransfrag> keeptrf;
	GVec<int> checktrf;
	for(int f=trflong.Count()-1;f>=0;f--) { // if this is a guide it should be reflected in the prediction downstream
		path.Clear();
		 int t=trflong[f];
		 if(t<0) GError("Stored long transcript is negative!\n");
		 pathpat=transfrag[t]->pattern;
		 minp=transfrag[t]->nodes[0];
		 maxp=transfrag[t]->nodes.Last();

		 int *pos=gpos[edge(0,minp,gno)];
		 if(pos) pathpat[*pos]=1;

		 pos=gpos[edge(maxp,gno-1,gno)];
		 if(pos) pathpat[*pos]=1;

		 maxi=minp;
		 path.Add(maxi);
		 pathpat[maxi]=1;

		 istranscript.reset();

		 float flux=0;
		 GVec<float> nodeflux;

		 /*
	 	 { // DEBUG ONLY
	 	 fprintf(stderr,"\n\n***Start get_trf_long_mix with maxi=%d minp=%d maxp=%d guide=%d and transcript:",maxi,minp,maxp,transfrag[t]->guide);
	 	 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
	 	 fprintf(stderr," pathpat=");
	 	 //printBitVec(pathpat);
		 fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf_unispg memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 	 }
	 	 */


		 bool tocheck=true;
		 if(back_to_source_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
			 path.cAdd(0);
			 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			 if(fwd_to_sink_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {

				 flux=long_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);

				 /*
				 { // DEBUG ONLY
					 //printTime(stderr);
					 fprintf(stderr,"flux=%g Path:",flux);
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
					 fprintf(stderr,"\n");
					 fprintf(stderr,"Nodecapacities:");
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %f",nodeflux[i]);
					 fprintf(stderr,"***\n");
				 }
				 */

				 if(flux) { // these are not valid paths in the graph


					 tocheck=false;

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 int startnode=1;
					 int lastnode=path.Count()-2;

					 uint startpoint=no2gnode[path[startnode]]->end;
					 uint endpoint=no2gnode[path[lastnode]]->start;
					 CGraphnode *jnode=no2gnode[path[startnode]];
					 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through startnode
						 int t=jnode->trf[i];
						 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
							 if(transfrag[t]->nodes[0]==path[startnode]) {
								 if(transfrag[t]->longstart) {
									 if(transfrag[t]->longstart<startpoint) { startpoint=transfrag[t]->longstart;}
								 }
							 }
						 }
					 }
					 jnode=no2gnode[path[lastnode]];
					 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through lastnode
						 int t=jnode->trf[i];
						 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
							 if(transfrag[t]->nodes.Last()==path[lastnode]) {
								 if(transfrag[t]->longend) {
									 if(transfrag[t]->longend>endpoint) { endpoint=transfrag[t]->longend;}
								 }
							 }
						 }
					 }
					 if(startpoint==no2gnode[path[startnode]]->end) startpoint=no2gnode[path[1]]->start;
					 if(endpoint==no2gnode[path[lastnode]]->start) endpoint=no2gnode[path[lastnode]]->end;


					 while(j<=lastnode) {
						 int nodestart=no2gnode[path[j]]->start;
						 int nodeend=no2gnode[path[j]]->end;
						 len+=nodeend-nodestart+1;
						 nodecov[path[j]]-=nodeflux[j];
						 float ecov=nodeflux[j]*noderate[path[j]];
						 float excov=ecov;
						 /*if(mixedMode) {
							 no2gnode[path[j]]->cov-=ecov;
							 if(no2gnode[path[j]]->cov<epsilon) no2gnode[path[j]]->cov=0;
						 }*/
						 //float excov=nodeflux[j]*noderate[path[j]];
						 while(j+1<=lastnode && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 //while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 nodeend=no2gnode[path[j]]->end;
							 ecov=nodeflux[j]*noderate[path[j]];
							 nodecov[path[j]]-=nodeflux[j];
							 excov+=ecov;
							 len+=nodeend-no2gnode[path[j]]->start+1;
						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 cov+=excov;
						 exoncov.Add(excov);
						 j++;
					 }
					 if(transfrag[t]->nodes.Count()==1) transfrag[t]->abundance=0;
					 //fprintf(stderr,"Store prediction %d  with abundance=%f len=%d\n",pred.Count(),cov/len,len);
					 //GffObj *g=NULL;

					 if(len>=mintranscriptlen) {
						 if(first) { geneno++; first=false;}
						 /*fprintf(stderr,"1 Store prediction %d  with abundance=%f totalabundance=%f len=%d startpoint=%d endpoint=%d and exons:",pred.Count(),cov/len,cov,len,startpoint,endpoint);
						 for(int i=0;i<exons.Count();i++) fprintf(stderr," %d-%d",exons[i].start,exons[i].end);
						 fprintf(stderr,"\n");*/
						 GffObj *g=NULL;
						 if(transfrag[t]->guide) {
							 g=guides[int(transfrag[t]->guide-1)];
							 if (g && g->uptr) {
								 RC_TData &td = *(RC_TData*) (g->uptr);
								 td.in_bundle=3;
								 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
							 }
						 }

						 CPrediction *p=new CPrediction(geneno, g,startpoint , endpoint, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->mergename='.'; // I should not delete this prediction
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 //fprintf(stderr,"Added prediction=%d with totalcov=%.1f\n",pred.Count()-1,pred.Last()->cov);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
			 }
		 }

		 if(tocheck)  { // try to see if you can rescue transfrag
			if(!guided || transfrag[t]->guide || (no2gnode[transfrag[t]->nodes[0]]->parent[0]==0 && 
				   no2gnode[transfrag[t]->nodes.Last()]->child.Last()==gno-1) )
				// only accept long transfrags that are linked to source and sink
			 checktrf.Add(t);
		 }
	 }

	//keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete

	 for(int c=0;c<checktrf.Count();c++) if(transfrag[checktrf[c]]->guide || transfrag[checktrf[c]]->abundance>=readthr) { // only in this case it is worth considering it as a potential prediction
		 int t=checktrf[c];

		 /*fprintf(stderr,"checktrf[%d]=%d with abundance=%f with start=%d end=%d nodes:",c,t,transfrag[t]->abundance,transfrag[t]->longstart,transfrag[t]->longend);
		 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
		 fprintf(stderr,"\n");*/

		 GVec<int> tmatch;
		 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch);
		 if(abundancesum>0) {
			 if(!transfrag[t]->shortread && transfrag[t]->nodes.Count()>1) {
				 for(int j=0;j<tmatch.Count();j++) { // found good match(es) but transcript has more than one node for this to make sense
					 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum;
					 int p=0;
					 int i=0;
					 int np=npred+tmatch[j]; // how do I know that keeptrf lead to a prediction -> because keeptrf represent the predictions that were added
					 //fprintf(stderr,"Add %.1f to prediction %d with cov=%.1f\n",abundprop,np,pred[np]->cov);
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 //fprintf(stderr,"newnodecov=%.1f\n",newnodecov);
								 if(newnodecov<0) {
									 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
									 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
									 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;

								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
								 //fprintf(stderr,"...add cov=%1.f (totalcov=%.1f) to exon[%d] and pred[%d]->cov=%f\n",addcov/no2gnode[transfrag[t]->nodes[i]]->len(),addcov,p,np,pred[np]->cov);
							 }
							 i++;
						 }
					 }
				 }
				 transfrag[t]->abundance=0;
			 }
		 }
		 else {
			 if(!eonly || transfrag[t]->guide) { // store it as an independent prediction
				 pathpat=transfrag[t]->pattern; // not used right now but maybe in the future?
				 path.Clear();
				 path.cAdd(0);
				 path.Add(transfrag[t]->nodes[0]);
				 for(int j=1;j<transfrag[t]->nodes.Count();j++) {
					 if(transfrag[t]->nodes[j]!=1+transfrag[t]->nodes[j-1] ||
							 no2gnode[transfrag[t]->nodes[j]]->start-1!=no2gnode[transfrag[t]->nodes[j-1]]->end) {
						 // check if transfrag t1 is incomplete between node[n-1] and node [n]
						 int *pos=gpos[edge(transfrag[t]->nodes[j-1],transfrag[t]->nodes[j],gno)];
						 if(!pos || !transfrag[t]->pattern[*pos]) { // incomplete transfrag
							 break;
						 }
						 if(pos) pathpat[*pos]=1;
						 path.Add(transfrag[t]->nodes[j]);
					 }
					 else path.Add(transfrag[t]->nodes[j]);
				 }
				 if(path.Last()==transfrag[t]->nodes.Last()) { // this transfrag is complete, might be worth rescuing

					 uint startpoint=no2gnode[path[1]]->end;
					 uint endpoint=no2gnode[path.Last()]->start;
					 if(transfrag[t]->longstart) startpoint=transfrag[t]->longstart;
					 if(transfrag[t]->longend) endpoint=transfrag[t]->longend;

					 int sink=gno-1;
					 path.Add(sink);

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 while(j<path.Count()-1) {
						 int nodestart=no2gnode[path[j]]->start;
						 int nodeend=no2gnode[path[j]]->end;
						 // nodecov[path[j]]-=transfrag[t]->abundance; // do not need this here anymore
						 len+=nodeend-nodestart+1;

						 float excov=0;
						 if(nodecov[path[j]]>epsilon) {
							 excov=transfrag[t]->abundance*no2gnode[path[j]]->len();
							 float newnodecov=nodecov[path[j]]-excov/noderate[path[j]];
							 if(newnodecov<0) {
								 excov=nodecov[path[j]]*noderate[path[j]];
								 newnodecov=0;
							 }
							 nodecov[path[j]]=newnodecov;
						 }

						 while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 len+=no2gnode[path[j]]->len();
							 nodeend=no2gnode[path[j]]->end;

							 float addcov=0;
							 if(nodecov[path[j]]>epsilon) {
								 addcov=transfrag[t]->abundance*no2gnode[path[j]]->len();
								 float newnodecov=nodecov[path[j]]-addcov/noderate[path[j]];
								 if(newnodecov<0) {
									 addcov=nodecov[path[j]]*noderate[path[j]];
									 newnodecov=0;
								 }
								 nodecov[path[j]]=newnodecov;
							 }
							 excov+=addcov;
							 //excov+=transfrag[t]->abundance*noderate[path[j]];
						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 cov+=excov;
						 exoncov.Add(excov);
						 j++;
					 }
					 //GffObj *g=NULL;
					 if(len>=mintranscriptlen) {
						 if(first) { geneno++; first=false;}
						 GffObj *g=NULL;
						 if(transfrag[t]->guide) {
							 g=guides[int(transfrag[t]->guide-1)];
							 if (g && g->uptr) {
								 RC_TData &td = *(RC_TData*) (g->uptr);
								 td.in_bundle=3;
								 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
							 }
						 }
						 //fprintf(stderr,"2 Store prediction %d:%d-%d  with len=%d and abundance=%f startpoint=%d endpoint=%d\n",pred.Count(),exons[0].start ,exons.Last().end,len,cov/len,startpoint,endpoint);
						 CPrediction *p=new CPrediction(geneno, g,startpoint , endpoint, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
			 }
			 transfrag[t]->abundance=0;
		 }
	 }


	 if(pred.Count()>npred) {
		 for(int t=0;t<transfrag.Count();t++) if(transfrag[t]->longread) {
			 if(!transfrag[t]->shortread && transfrag[t]->nodes.Count()>1 && transfrag[t]->abundance>epsilon && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1 ) {
				 /*fprintf(stderr,"Consider transfrag[%d]->abundance=%f with start=%d end=%d nodes:",t,transfrag[t]->abundance,transfrag[t]->longstart,transfrag[t]->longend);
				 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
				 fprintf(stderr,"\n");*/

				 GVec<int> tmatch;
				 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch); // abundancesum is the sum of all matching transcripts
				 if(abundancesum>0) for(int j=0;j<tmatch.Count();j++){ // found good match(es)
					 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum; // proportion of transcript that will be allocated to this matching one
					 int p=0;
					 int i=0;
					 int np=npred+tmatch[j];
					 //fprintf(stderr,"Add %.1f to prediction %d with cov=%.1f",abundprop,np,pred[np]->cov);
					 //for(int i=0;i<keeptrf[tmatch[j]].nodes.Count();i++) fprintf(stderr," %d",keeptrf[tmatch[j]].nodes[i]);
					 //fprintf(stderr,"\n");
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 //fprintf(stderr,"newnodecov=%.1f\n",newnodecov);
								 if(newnodecov<0) {
									 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
									 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
									 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;
								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
								 //fprintf(stderr,"...add cov=%1.f (totalcov=%.1f) to exon[%d] and pred[%d]->cov=%f\n",addcov/no2gnode[transfrag[t]->nodes[i]]->len(),addcov,p,np,pred[np]->cov);
							 }
							 i++;
						 }
					 }
				 }
			 }
			 if(!transfrag[t]->nodes[0] || transfrag[t]->nodes.Last()==gno-1) // deletes starting transfrangs
				 transfrag[t]->abundance=0;
			 else transfrag[t]->abundance=transfrag[t]->usepath; // this restores tranfrag[t] abundance to what it was before so when I do the short read flux I can utilizes fully

			 transfrag[t]->usepath=-1;
		 }

		 int p=npred;
		 while(p<pred.Count()) {
			 if(pred[p]->cov) {
				 //fprintf(stderr,"Pred[%d] has coverage=%.1f (totalcov=%.1f)\n",p,pred[p]->cov,pred[p]->cov/abs(pred[p]->tlen));
				 pred[p]->cov/=abs(pred[p]->tlen);
				 for(int i=0;i<pred[p]->exons.Count();i++)
					 pred[p]->exoncov[i]/=pred[p]->exons[i].len();
				 // adjust start/endpoints
				 if(pred[p]->start!=pred[p]->exons[0].start) {
					 pred[p]->tlen+=pred[p]->start-pred[p]->exons[0].start;
					 pred[p]->exons[0].start=pred[p]->start;
				 }
				 if(pred[p]->end!=pred[p]->exons.Last().end) {
					 pred[p]->tlen+=pred[p]->exons.Last().end-pred[p]->end;
					 pred[p]->exons.Last().end=pred[p]->end;
				 }
				 p++;
			 }
			 else if(!eonly) { // || !pred[p]->t_eq) {
				 //fprintf(stderr,"delete prediction %d\n",p);
				 pred.Delete(p); // I delete all predictions that have 0 coverage unless it's eonly mode
			 }
			 else p++;
		 }
	 }

	 int nkeep=keeptrf.Count();
	 /*for(int i=nkeep-1;i>nkept-1;i--) { // delete all kept transfrags we are not confident in
			 keeptrf.Delete(i);
	 }*/
	 for(int i=0;i<nkeep;i++) keeptrf[i].weak=i+npred; // I need to remember the prediction it represents
	 keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph comes first, then the one with most nodes, then the one more complete

	 float flux=0;
	 GVec<float> nodeflux;
	 path.Clear();


	 for(int i=0;i<nkeep;i++) { // compute flux from short read data here
		 istranscript.reset();
		 nodeflux.Clear();
		 bool full=true;
		 // make sure source and sink are present and they are connected in the pattern of the prediction
		 if(keeptrf[i].nodes[0]) {
			 keeptrf[i].nodes.Insert(0,0);
			 keeptrf[i].pattern[0]=1;
		 }
		 int sink=gno-1;
		 if(keeptrf[i].nodes.Last()!=gno-1) {
			 keeptrf[i].nodes.Add(sink);
			 keeptrf[i].pattern[sink]=1;
		 }
		 int key=edge(0,keeptrf[i].nodes[1],gno);
		 int *pos=NULL;
		 if(key<(int)gpos.Count()) pos=gpos[key];
		 if(pos!=NULL) keeptrf[i].pattern[*pos]=1;
		 pos=NULL;
		 int n=keeptrf[i].nodes.Count();
		 key=edge(keeptrf[i].nodes[n-2],sink,gno);
		 if(pos!=NULL) keeptrf[i].pattern[*pos]=1;

		 flux=push_max_flow(gno,keeptrf[i].nodes,istranscript,transfrag,no2gnode,nodeflux,keeptrf[i].pattern,gpos,full);

		 /*
		 { // DEBUG ONLY
			 //printTime(stderr);
			 fprintf(stderr,"flux=%g Path:",flux);
			 for(int j=0;j<keeptrf[i].nodes.Count();j++) fprintf(stderr," %d",keeptrf[i].nodes[j]);
			 fprintf(stderr,"***\n");
		 }
		 */

		 if(flux>epsilon) {

			 bool included=true;
			 float cov=store_transcript(pred,keeptrf[i].nodes,nodeflux,nodecovall,no2gnode,geneno,first,strand,gno,gpos,included,prevpath);
			 if(cov) {

				 //fprintf(stderr,"Coverage of long pred[%d]=%f vs. short pred=%f\n",keeptrf[i].weak,pred[keeptrf[i].weak]->cov,pred.Last()->cov);

				 if(pred.Last()->cov>pred[keeptrf[i].weak]->cov) { // new prediction is better than the previous one -- shouldn't I add this to previous prediction?
					 int p=keeptrf[i].weak;
					 pred[p]->cov=pred.Last()->cov;
					 for(int k=0;k<pred[p]->exons.Count();k++) {
						 pred[p]->exoncov[k]=pred.Last()->exoncov[k];
					 }
					 /*pred[p]->start=pred.Last()->start;
					 pred[p]->end=pred.Last()->end;
					 pred[p]->exons[0].start=pred.Last()->exons[0].start;
					 pred[p]->exons.Last().end=pred.Last()->exons.Last().end;
					 pred[p]->tlen=pred.Last()->tlen;*/
					 //pred[p]->tlen=abs(pred[p]->tlen);
				 }
				 delete pred.Pop(); //prevent memory leak, popped element is otherwise "forgotten"
			 }
		 }

	 }

}

void get_trf_long_unispg(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,int strand,GList<CPrediction>& pred,GVec<int>& trflong,BundleData *bdata) {

	GPVec<GffObj>& guides = bdata->keepguides;

	GVec<float> nodecov; // the coverage of all transfrags entering a node
	GVec<float> noderate;
	for(int i=0;i<gno;i++) {
		CGraphnode *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0); // this are all transfrags that link nodes together
		float rate=1;
		if(i && i<gno-1) {
			for(int j=0;j<inode->trf.Count();j++) { // for all transfrags going through node
				int t=inode->trf[j];
				//if(transfrag[t]->nodes[0]<i) { // entering transfrags:
				if(transfrag[t]->nodes.Last()>i) { // exiting transfrags: this is more consistent with the nodeflux computation
					nodecov[i]+=transfrag[t]->abundance;
				}
			}

			if(nodecov[i]) rate=nodecov[i];
			if(rate<=0) rate=1; // this shouldn't happen
			//fprintf(stderr,"rate=%f\n",rate);
			rate=inode->cov/rate;
		}
		noderate.Add(rate);
		//fprintf(stderr,"Node[%d]:%d-%d no2gnode->cov=%f nodecov=%f noderate=%f\n",i,inode->start,inode->end,inode->cov,nodecov[i],noderate[i]);
	}

	GBitVec istranscript(transfrag.Count());

	 GVec<int> path;
	 GBitVec pathpat(gno+edgeno);
	 int minp;
	 int maxp;
	 int maxi;

	 char sign='-';
	 if(strand) { sign='+';}
	 int npred=pred.Count();

	 GVec<CTransfrag> keeptrf;
	 GVec<int> checktrf;
	 for(int f=trflong.Count()-1;f>=0;f--) { // if this is a guide it should be reflected in the prediction downstream
		 path.Clear();
		 int t=trflong[f];
		 if(t<0) GError("Stored long transcript is negative!\n");
		 pathpat=transfrag[t]->pattern;
		 minp=transfrag[t]->nodes[0];
		 maxp=transfrag[t]->nodes.Last();

		 //if(no2gnode[transfrag[t]->nodes[0]]->hardstart) {
			 int *pos=gpos[edge(0,minp,gno)];
			 if(pos) pathpat[*pos]=1;
			 //minp=0;
		 //}
		 //if(no2gnode[transfrag[t]->nodes.Last()]->hardend) {
			 pos=gpos[edge(maxp,gno-1,gno)];
			 if(pos) pathpat[*pos]=1;
			 //maxp=gno-1;
		 //}

		 maxi=minp;
		 path.Add(maxi);
		 pathpat[maxi]=1;

		 istranscript.reset();

		 float flux=0;
		 //float fragno=0;
		 GVec<float> nodeflux;

		 /*
	 	 { // DEBUG ONLY
	 	 fprintf(stderr,"\n\n***Start get_trf_long with maxi=%d minp=%d maxp=%d guide=%d and transcript:",maxi,minp,maxp,transfrag[t]->guide);
	 	 for(int i=0;i<transfrag[t]->nodes.Count();i++) {
	 		 if(i) {
	 			 pos=gpos[edge(transfrag[t]->nodes[i-1],transfrag[t]->nodes[i],gno)];
	 			 if(pos && pathpat[*pos])
	 				 fprintf(stderr,"-");
	 		 }
	 		 fprintf(stderr," %d",transfrag[t]->nodes[i]);

	 	 }
	 	 fprintf(stderr," pathpat=");
	 	 //printBitVec(pathpat);
		 fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf_unispg memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 	 }
	 	 */

		 bool tocheck=true;
		 if(back_to_source_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
			 path.cAdd(0);
			 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			 if(fwd_to_sink_fast_long(maxi,path,minp,maxp,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {

				 flux=long_max_flow(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat);

				 /*
				 { // DEBUG ONLY
					 //printTime(stderr);
					 fprintf(stderr,"flux=%g Path:",flux);
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %d",path[i]);
					 fprintf(stderr,"\n");
					 fprintf(stderr,"Nodecapacities:");
					 for(int i=0;i<path.Count();i++) fprintf(stderr," %f",nodeflux[i]);
					 fprintf(stderr,"***\n");
				 }
				 */

				 if(flux) { // these are not valid paths in the graph
					 tocheck=false;

					 GVec<GSeg> exons;
					 GVec<float> exoncov;
					 int j=1;
					 int len=0;
					 float cov=0;
					 //int startnode=j;
					 int lastnode=path.Count()-2;
					 //uint startpoint=no2gnode[path[lastnode]]->end;
					 //uint endpoint=no2gnode[path[1]]->start;

					 /* if(mixedMode) { // establish start/end point of path
						 while(j<path.Count()-1) {
							 CGraphnode *jnode=no2gnode[path[j]];
							 for(int i=0;i<jnode->trf.Count();i++) { // for all transfrags going through node
								 int t=jnode->trf[i];
								 if(istranscript[t] && transfrag[t]->longread && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) {
									 if(transfrag[t]->nodes[0]==path[j]) {
										 if(transfrag[t]->longstart) {
											 if(transfrag[t]->longstart<startpoint) { startpoint=transfrag[t]->longstart; startnode=j;}
										 }
										 else if(jnode->start<startpoint) { startpoint=jnode->start; startnode=j;}
									 }
									 if(transfrag[t]->nodes.Last()==path[j]) {
										 if(transfrag[t]->longend) {
											 if(transfrag[t]->longend>endpoint) { endpoint=transfrag[t]->longend; lastnode=j;}
										 }
										 else if(jnode->end>endpoint) {endpoint=jnode->end; lastnode=j;}
									 }
								 }
							 }
							 j++;
						 }
						 if(startpoint>endpoint) {
							 j=1;
							 lastnode=path.Count()-2;
							 startpoint=no2gnode[path[1]]->start;
							 endpoint=no2gnode[path[lastnode]]->end;
						 }
						 else j=startnode;
						 //fprintf(stderr,"startnode=%d lastnode=%d startpoint=%d endpoint=%d\n",path[startnode],path[lastnode],startpoint,endpoint);
					 }*/

					 while(j<=lastnode) {
					 //while(j<path.Count()-1) {
						 int nodestart=no2gnode[path[j]]->start;
						 int nodeend=no2gnode[path[j]]->end;
						 nodecov[path[j]]-=nodeflux[j];
						 len+=nodeend-nodestart+1;
						 float ecov=nodeflux[j]*noderate[path[j]];
						 float excov=ecov;
						 //fprintf(stderr,"excov+=%f * %f = %f\n",nodeflux[j],noderate[path[j]],excov);
						 /*if(mixedMode) {
							 no2gnode[path[j]]->cov-=ecov;
							 if(no2gnode[path[j]]->cov<epsilon) no2gnode[path[j]]->cov=0;
						 }*/
						 //float excov=nodeflux[j]*noderate[path[j]];
						 while(j+1<=lastnode && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 //while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
							 j++;
							 nodecov[path[j]]-=nodeflux[j];
							 nodeend=no2gnode[path[j]]->end;
							 ecov=nodeflux[j]*noderate[path[j]];
							 excov+=ecov;
							 //fprintf(stderr,"excov+=%f * %f = %f\n",nodeflux[j],noderate[path[j]],excov);
							 len+=nodeend-no2gnode[path[j]]->start+1;

						 }
						 GSeg exon(nodestart,nodeend);
						 exons.Add(exon);
						 //fprintf(stderr,"excov=%f\n",excov/(exon.end-exon.start+1));
						 cov+=excov;
						 //fprintf(stderr,"cov+=%f=%f\n",excov,cov);
						 exoncov.Add(excov);
						 j++;
					 }
					 if(transfrag[t]->nodes.Count()==1) transfrag[t]->abundance=0;
					 //fprintf(stderr,"Store prediction %d  with abundance=%f len=%d\n",pred.Count(),cov/len,len);
					 GffObj *g=NULL;
					 if(transfrag[t]->guide) {
						 g=guides[int(transfrag[t]->guide-1)];
						 if (g && g->uptr) {
							 RC_TData &td = *(RC_TData*) (g->uptr);
							 td.in_bundle=3;
							 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
						 }
					 }
					 if(!eonly || g) {
						 /*fprintf(stderr,"1 Store prediction %d  with abundance=%f len=%d and exons:",pred.Count(),cov/len,len);
						 for(int i=0;i<exons.Count();i++) fprintf(stderr," %d-%d",exons[i].start,exons[i].end);
						 fprintf(stderr,"\n");*/
						 CPrediction *p=new CPrediction(geneno, g,exons[0].start , exons.Last().end, cov, sign, len);
						 p->exons=exons;
						 p->exoncov=exoncov;
						 p->mergename='.'; // I should not delete this prediction
						 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
						 pred.Add(p);

						 CTransfrag u(path,pathpat,cov/len);
						 keeptrf.Add(u);
					 }
				 }
				 else if(transfrag[t]->guide) {
					 checktrf.Add(t);
				 }
			 }
		 }

		 if(tocheck)  { // try to see if you can rescue transfrag -> they are stored from more abundant to least -> not if using mixedMode
			 checktrf.Add(t);
		 }
	 }

	 //keeptrf.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete

	 for(int c=0;c<checktrf.Count();c++) if(transfrag[checktrf[c]]->guide || transfrag[checktrf[c]]->abundance>=readthr) { // only in this case it is worth considering it as a potential prediction
		 int t=checktrf[c];

		 //fprintf(stderr,"checktrf[%d]=%d with abundance=%f\n",c,t,transfrag[t]->abundance);

		 GVec<int> tmatch;
		 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch);
		 if(abundancesum>0) for(int j=0;j<tmatch.Count();j++){ // found good match(es)
			 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum;
			 int p=0;
			 int i=0;
			 int np=npred+tmatch[j]; // how do I know that keeptrf lead to a prediction -> because keeptrf represent the predictions that were added
			 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
				 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
				 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
				 else { // the two intersect (I can only have the full node included in exon)
					 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
						 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
						 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
						 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
						 if(newnodecov<0) {
							 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
							 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
							 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
						 }
						 nodecov[transfrag[t]->nodes[i]]=newnodecov;

						 pred[np]->exoncov[p]+=addcov;
						 pred[np]->cov+=addcov;
					 }
					 i++;
				 }
			 }
		 }
		 else if(!eonly || transfrag[t]->guide) { // store it as an independent prediction
			 pathpat=transfrag[t]->pattern; // not used right now but maybe in the future?
			 path.Clear();
			 path.cAdd(0);
			 path.Add(transfrag[t]->nodes[0]);
			 for(int j=1;j<transfrag[t]->nodes.Count();j++) {
				 if(transfrag[t]->nodes[j]!=1+transfrag[t]->nodes[j-1] ||
						 no2gnode[transfrag[t]->nodes[j]]->start-1!=no2gnode[transfrag[t]->nodes[j-1]]->end) {
					 // check if transfrag t1 is incomplete between node[n-1] and node [n]
					 int *pos=gpos[edge(transfrag[t]->nodes[j-1],transfrag[t]->nodes[j],gno)];
					 if(!pos || !transfrag[t]->pattern[*pos]) { // incomplete transfrag
						 break;
					 }
					 if(pos) pathpat[*pos]=1;
					 path.Add(transfrag[t]->nodes[j]);
				 }
				 else path.Add(transfrag[t]->nodes[j]);
			 }
			 if(path.Last()==transfrag[t]->nodes.Last()) { // this transfrag is complete, might be worth rescuing
				 int sink=gno-1;
				 path.Add(sink);

				 GVec<GSeg> exons;
				 GVec<float> exoncov;
				 int j=1;
				 int len=0;
				 float cov=0;
				 while(j<path.Count()-1) {
					 int nodestart=no2gnode[path[j]]->start;
					 int nodeend=no2gnode[path[j]]->end;
					 // nodecov[path[j]]-=transfrag[t]->abundance; // do not need this here anymore
					 len+=nodeend-nodestart+1;

					 float excov=0;
					 if(nodecov[path[j]]>epsilon) {
						 excov=transfrag[t]->abundance*no2gnode[path[j]]->len();
						 float newnodecov=nodecov[path[j]]-excov/noderate[path[j]];
						 if(newnodecov<0) {
							 excov=nodecov[path[j]]*noderate[path[j]];
							 newnodecov=0;
						 }
						 nodecov[path[j]]=newnodecov;
					 }

					 while(j+1<path.Count()-1 && no2gnode[path[j]]->end+1==no2gnode[path[j+1]]->start) {
						 j++;
						 len+=no2gnode[path[j]]->len();
						 nodeend=no2gnode[path[j]]->end;

						 float addcov=0;
						 if(nodecov[path[j]]>epsilon) {
							 addcov=transfrag[t]->abundance*no2gnode[path[j]]->len();
							 float newnodecov=nodecov[path[j]]-addcov/noderate[path[j]];
							 if(newnodecov<0) {
								 addcov=nodecov[path[j]]*noderate[path[j]];
								 newnodecov=0;
							 }
							 nodecov[path[j]]=newnodecov;
						 }
						 excov+=addcov;
						 //excov+=transfrag[t]->abundance*noderate[path[j]];
					 }
					 GSeg exon(nodestart,nodeend);
					 exons.Add(exon);
					 cov+=excov;
					 exoncov.Add(excov);
					 j++;
				 }
				 GffObj *g=NULL;
				 if(transfrag[t]->guide) {
					 g=guides[int(transfrag[t]->guide-1)];
					 if (g && g->uptr) {
						 RC_TData &td = *(RC_TData*) (g->uptr);
						 td.in_bundle=3;
						 //fprintf(stderr,"sg guide %s is stored\n",g->getID());
					 }
				 }
				 //fprintf(stderr,"2 Store prediction %d:%d-%d  with len=%d and abundance=%f\n",pred.Count(),exons[0].start ,exons.Last().end,len,cov/len);
				 CPrediction *p=new CPrediction(geneno, g,exons[0].start , exons.Last().end, cov, sign, len);
				 p->exons=exons;
				 p->exoncov=exoncov;
				 p->tlen=-p->tlen; // negative transcript length signifies assembly is from a long read
				 pred.Add(p);

				 CTransfrag u(path,pathpat,cov/len);
				 keeptrf.Add(u);
			 }
		 }
		 transfrag[t]->abundance=0;
	 }

	 if(pred.Count()>npred) {
		 /*if(mixedMode) {
			 for(int t=0;t<transfrag.Count();t++) {
			   if(transfrag[t]->longread && (!transfrag[t]->nodes[0] || transfrag[t]->nodes.Last()==gno-1) )
			     transfrag[t]->abundance=0;
			   else
			     transfrag[t]->abundance=transfrag[t]->usepath;
			   transfrag[t]->usepath=-1;
			 }
		 }
		 else */for(int t=0;t<transfrag.Count();t++) if(transfrag[t]->longread) { // longreads mode --> tries to add all transfrags to predictions
			 if(transfrag[t]->abundance>epsilon && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1 ) {
			   /*fprintf(stderr,"Consider transfrag[%d]->abundance=%f with start=%d end=%d nodes:",t,transfrag[t]->abundance,transfrag[t]->longstart,transfrag[t]->longend);
				 for(int i=0;i<transfrag[t]->nodes.Count();i++) fprintf(stderr," %d",transfrag[t]->nodes[i]);
				 fprintf(stderr,"\n");*/

				 GVec<int> tmatch;
				 float abundancesum=best_trf_match(transfrag[t],keeptrf,no2gnode,gno,tmatch);
				 if(abundancesum>0) for(int j=0;j<tmatch.Count();j++){ // found good match(es)
					 float abundprop=transfrag[t]->abundance*keeptrf[tmatch[j]].abundance/abundancesum;
					 //fprintf(stderr,"Added %.1f to prediction:",abundprop);
					 //for(int i=0;i<keeptrf[tmatch[j]].nodes.Count();i++) fprintf(stderr," %d",keeptrf[tmatch[j]].nodes[i]);
					 int p=0;
					 int i=0;
					 int np=npred+tmatch[j];
					 while(i<transfrag[t]->nodes.Count() && p<pred[np]->exons.Count()) {
						 if(no2gnode[transfrag[t]->nodes[i]]->end<pred[np]->exons[p].start) i++;
						 else if(pred[np]->exons[p].end<no2gnode[transfrag[t]->nodes[i]]->start) p++;
						 else { // the two intersect (I can only have the full node included in exon)
							 if(nodecov[transfrag[t]->nodes[i]]>epsilon){
								 //float addcov=transfrag[t]->abundance*noderate[transfrag[t]->nodes[i]];
								 float addcov=abundprop*no2gnode[transfrag[t]->nodes[i]]->len();
								 float newnodecov=nodecov[transfrag[t]->nodes[i]]-addcov/noderate[transfrag[t]->nodes[i]];
								 if(newnodecov<0) {
									 newnodecov=nodecov[transfrag[t]->nodes[i]]*keeptrf[tmatch[j]].abundance/abundancesum;
									 addcov=newnodecov*noderate[transfrag[t]->nodes[i]];
									 newnodecov=nodecov[transfrag[t]->nodes[i]]-newnodecov;
								 }
								 nodecov[transfrag[t]->nodes[i]]=newnodecov;
								 //fprintf(stderr," newnodecov=%f\n",newnodecov);
								 pred[np]->exoncov[p]+=addcov;
								 pred[np]->cov+=addcov;
							 }
							 i++;
						 }
					 }
					 //fprintf(stderr," new abundance=%f\n",pred[np]->cov/pred[np]->tlen);
				 }
			 }
			 transfrag[t]->abundance=0; // delete abundance in order not to use it in short reads
		 }

		 int p=npred;
		 while(p<pred.Count()) {
			 if(pred[p]->cov) {
				 pred[p]->cov/=abs(pred[p]->tlen);
				 for(int i=0;i<pred[p]->exons.Count();i++)
					 pred[p]->exoncov[i]/=pred[p]->exons[i].len();
				 p++;
			 }
			 else if(!eonly) { // || !pred[p]->t_eq) {
				 //fprintf(stderr,"delete prediction %d\n",p);
				 pred.Delete(p); // I delete all predictions that have 0 coverage unless it's eonly mode
			 }
			 else p++;
		 }
	 }
}

