#include "findTranscripts_A.h"

// int find_transcripts(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand, BundleData* bdata) {

int find_transcripts_APPLY_UNISPG(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,int geneno,int strand, GVec<CGuide>& guidetrf,GPVec<GffObj>& guides,GVec<int>& guidepred,BundleData* bdata,GVec<int>& trflong) {
	GList<CPrediction>& pred = bdata->pred;
	fprintf(stderr, "find_transcripts_APPLY_UNISPG\n");

	// process in and out coverages for each node
	int maxi=0; // node with maximum coverage
	GVec<float> nodecov; // node coverages

	for(int i=0;i<gno;i++) {
		fprintf(stderr, "i: %d\n", i);
		CGraphnodeUnispg *inode=no2gnode[i]; // this is here only because of the DEBUG option below
		nodecov.cAdd(0.0);

		if(i) { // for all nodes but the source

		    if(i<gno-1 && inode->len()) nodecov[i]=(inode->cov_s.Get(0))/inode->len(); // sink also has 0 coverage

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

		    if(abundin) inode->rate_s[0]=abundout/abundin;
		    if(abundout) inode->capacity_s[0]=abundout+abundthrough; // node capacity tells me how much of that node coverage I can use given how many transfrags leave the node
		    else inode->capacity_s[0]=abundin+abundthrough;

			// fprintf(stderr, "abundout: %f;  abundin: %f\n", abundout, abundin);
			// fprintf(stderr, "inode->capacity_s[0]: %f\n", inode->capacity_s.Get(0));
			fprintf(stderr, "inode->rate_s[0]: %f;  inode->capacity_s[0]: %f;  inode->capacity_s[0]: %f\n", inode->rate_s[0], inode->capacity_s[0], inode->capacity_s[0]);
		} // end if i

		/*
		{ // DEBUG ONLY
			printTime(stderr);
			fprintf(stderr,"Node %d: cov=%f capacity=%f rate=%f ",i,inode->cov/(inode->end-inode->start+1),inode->capacity,inode->rate);
			fprintf(stderr,"trf=");
			for(int t=0;t<inode->trf.Count();t++) fprintf(stderr," %d(%f)",inode->trf[t],transfrag[inode->trf[t]]->abundance);
			fprintf(stderr," maxi=%d maxcov=%f\n",maxi,nodecov[maxi]);
		}
		// */

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

	// process guides first
	//fprintf(stderr,"guidetrf.count=%d\n",guidetrf.Count());
	//if(guidetrf.Count()) maxi=guides_flow(gno,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat);

	bool first=true;

	fprintf(stderr,"guide count=%d\n",guidetrf.Count());

	if (eonly)
		guides_pushmaxflow_APPLY_UNISPG(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);
	else { //if(!eonly) {

		if(guidetrf.Count()) maxi=guides_pushmaxflow_APPLY_UNISPG(gno,edgeno,gpos,no2gnode,transfrag,guidetrf,geneno,strand,pred,nodecov,istranscript,pathpat,first,guides,guidepred,bdata);

		if(nodecov[maxi]>=1) { // sensitive mode only; otherwise >=readthr

			// 1:
			// parse_trf_weight_max_flow(gno,no2gnode,transfrag,geneno,strand,pred,nodecov,pathpat);
			// 2:
			fprintf(stderr, ">> In nodecov[maxi]: %f\n", nodecov[maxi]);
			GBitVec usednode(gno+edgeno);
			parse_trf_APPLY_UNISPG(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,usednode,0,pathpat);
		}
	}

    // /*
    { // DEBUG ONLY
    	for(int i=0;i<pred.Count();i++) {
    		if(pred[i]->t_eq) fprintf(stderr,"%s ",pred[i]->t_eq->getID());
    		fprintf(stderr,"pred[%d] (cov=%f,strand=%c):",i,pred[i]->cov,pred[i]->strand);
    		for(int j=0;j<pred[i]->exons.Count();j++) fprintf(stderr," %d-%d",pred[i]->exons[j].start,pred[i]->exons[j].end);
    		fprintf(stderr,"\n");
    	}
    }
    // */
	return(geneno);
}



int guides_pushmaxflow_APPLY_UNISPG(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata) {

	int maxi=1;
	int ng=guidetrf.Count();

	if(ng==1) { // if only one guide I do not need to do the 2 pass
		GVec<float> nodeflux;
		//float fragno=0;
		bool full=true;
		float flux= push_max_flow_APPLY_UNISPG(gno,guidetrf[0].trf->nodes,istranscript,transfrag,no2gnode,nodeflux,guidetrf[0].trf->pattern,gpos,full);
		istranscript.reset();

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"guide=%s flux[0]=%g\n",guides[guidetrf[0].g]->getID(),flux);
		}
		*/

		if(flux>epsilon) {
			bool include=true;
			if(guidepred[guidetrf[0].g]==-1) {

				store_transcript_APPLY_UNISPG(pred,guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[0].g]);
				//if(eonly) { // this is not correct because it might have been assigned before
					guidepred[guidetrf[0].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript_APPLY_UNISPG
					//fprintf(stderr,"guidepred[%d]=%d\n",guidetrf[0].g,guidepred[guidetrf[0].g]);
				//}
			}
			else {
				update_guide_pred_APPLY_UNISPG(pred,guidepred[guidetrf[0].g],guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
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
			guidetrf[g].trf->abundance=push_guide_maxflow_APPLY_UNISPG(gno,guidetrf[g].trf->nodes,istranscript,transfrag,no2gnode,guidetrf[g].trf->pattern);
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

			float flux=guidepushflow_APPLY_UNISPG(g,guidetrf,gno,istranscript,transfrag,no2gnode,nodeflux);

			istranscript.reset();

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"guide=%s flux[%d]=%f\n",guides[g]->getID(),g,flux);
			}
			*/

			bool include=true;
			if(flux>epsilon) {
				if(guidepred[guidetrf[g].g]==-1) {
					store_transcript_APPLY_UNISPG(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[g].g]);
					//if(eonly) {
						guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript_APPLY_UNISPG
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
					update_guide_pred_APPLY_UNISPG(pred,guidepred[guidetrf[g].g],guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
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
							store_transcript_APPLY_UNISPG(pred,guidetrf[g].trf->nodes,nodeflux,nodecov,no2gnode,geneno,first,s,gno,gpos,include,pathpat,false,bdata,guides[guidetrf[g].g]);
							//if(eonly) {
								guidepred[guidetrf[g].g]=pred.Count()-1; // NEED TO TEST: if this doesn't work for single genes I might want to recombine with the previous prediction in store_transcript_APPLY_UNISPG
								//fprintf(stderr,"2 guidepred[%d]=%d\n",guidetrf[g].g,guidepred[guidetrf[g].g]);
							//}
						}
						else {
							update_guide_pred_APPLY_UNISPG(pred,guidepred[guidetrf[g].g],guidetrf[0].trf->nodes,nodeflux,nodecov,no2gnode,gno,true);
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
				CTransfrag *trguide=find_guide_partial_pat_APPLY_UNISPG(guides[g],no2gnode,gno,edgeno,gpos,olen,olensum);
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
				if(guidepred[g]!=-1) update_guide_pred_APPLY_UNISPG(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
				else // new prediction for this guide
					guidepred[g]=store_guide_transcript_APPLY_UNISPG(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
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

				CGraphnodeUnispg *inode=no2gnode[n];
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
						if(guidepred[g]!=-1) update_guide_pred_APPLY_UNISPG(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,true);
						else // new prediction for this guide
							guidepred[g]=store_guide_transcript_APPLY_UNISPG(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],true);
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
									if(guidepred[g]!=-1) update_guide_pred_APPLY_UNISPG(pred,guidepred[g],path,nodeflux,nodecov,no2gnode,gno,false);
									else // new prediction for this guide
										guidepred[g]=store_guide_transcript_APPLY_UNISPG(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[g],false);
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
				if(guidepred[guidetrf[g].g]!=-1) update_guide_pred_APPLY_UNISPG(pred,guidepred[guidetrf[g].g],path,nodeflux,nodecov,no2gnode,gno,false);
				else // new prediction for this guide
					guidepred[guidetrf[g].g]=store_guide_transcript_APPLY_UNISPG(pred,path,nodeflux,nodecov,no2gnode,geneno,first,gno,guides[guidetrf[g].g],false);
			}
		}

		maxi=0;
	}

	return(maxi);
}



void parse_trf_APPLY_UNISPG(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,
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

	//  /*
	 { // DEBUG ONLY
	 	 fprintf(stderr,"\n\n***Start parse_trf with maxi=%d and cov=%f\n",maxi,nodecov[maxi]);
		 //fprintf(stderr,"Transcripts before path:");
		 //for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d",i);
		 //fprintf(stderr,"\n");

#ifdef GMEMTRACE
	 	 double vm,rsm;
	 	 get_mem_usage(vm, rsm);
	 	 GMessage("\t\tM(s):parse_trf memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
	 }
	//  */


	if(back_to_source_fast_APPLY_UNISPG(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
		fprintf(stderr, ">> back_to_source_fast_APPLY_UNISPG\n");
		 	 if(includesource) path.cAdd(0);
	 		 path.Reverse(); // back to source adds the nodes at the end to avoid pushing the list all the time

			if(fwd_to_sink_fast_APPLY_UNISPG(maxi,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos)) {
	 			 bool full=true;

				 flux=push_max_flow_APPLY_UNISPG(gno,path,istranscript,transfrag,no2gnode,nodeflux,pathpat,gpos,full);

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

	fprintf(stderr, ">> flux: %f\n", flux);
	 if(flux>epsilon) {
		 bool included=true;
		 float cov=store_transcript_APPLY_UNISPG(pred,path,nodeflux,nodecov,no2gnode,geneno,first,strand,gno,gpos,included,prevpath);

		//  /*
		 { // DEBUG ONLY
			 //fprintf(stderr,"Prevpath=");
			 //printBitVec(prevpath);
			 //fprintf(stderr,"\n");
		 	 fprintf(stderr,"cov=%f maxcov=%f\n",cov,maxcov);
		 }
		//  */

		 float frac=isofrac;
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
		 parse_trf_APPLY_UNISPG(maxi,gno,edgeno,gpos,no2gnode,transfrag,geneno,first,strand,pred,nodecov,istranscript,usednode,maxcov,prevpath);
	 }

}



// version of push_max_flow_APPLY_UNISPG where I weight the incomplete transfrags
float push_max_flow_APPLY_UNISPG(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,
		GVec<float>& nodeflux,GBitVec& pathpat, GIntHash<int> &gpos, bool &full) {

	int n=path.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
		nodeflux.cAdd(0.0);
	}
	GVec<float> capacityleft;   // how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n);
	capacityright.Resize(n);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n);
	sumright.Resize(n);

	//bool full=true;
	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"Start push max flow algorithm for path ");
		//printBitVec(pathpat);
		fprintf(stderr," :");
		for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		fprintf(stderr,"\n");
		//fprintf(stderr,"Used transcripts:");
		//for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		//fprintf(stderr,"\n");
	}
	*/

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance) {
				bool keeptr=false;
				if(istranscript[t]) keeptr=true;
				else if(!transfrag[t]->nodes[0]) {
					if(transfrag[t]->nodes.Last()==path[1]) keeptr=true;
				}
				else if(transfrag[t]->nodes.Last()==path.Last()) {
					if(transfrag[t]->nodes[0]==path[n-2]) keeptr=true;
				}
				else if(transfrag[t]->nodes[0]==path[i] && ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern)) { // only need to check transfrag the first time I encounter it
					keeptr=true;

					if(keeptr && !full) { // check if transcript fully supports path (full is false means I have not found any transcript to fully support path)
						full=true;
						int p=1;
						if(!transfrag[t]->longstart || !transfrag[t]->longend) full=false;
						if(full) while(path[p]<transfrag[t]->nodes[0]) {
							if(no2gnode[path[p]]->end+1<no2gnode[path[p+1]]->start) {
								full=false;
								break;
							}
							p++;
						}
						if(full) {
							p=path.Count()-2;
							while(path[p]>transfrag[t]->nodes.Last()) {
								if(no2gnode[path[p-1]]->end+1<no2gnode[path[p]]->start) {
									full=false;
									break;
								}
								p--;
							}
						}
						if(full) {
							for(p=2;p<path.Count()-2;p++) {
								if(!transfrag[t]->pattern[path[p]] && no2gnode[path[p]]->len()>longintronanchor) {
									full=false;
									break;
								}
							}
						}
					}
				}

				if(keeptr) { // transcript on path
					istranscript[t]=1;

					//fprintf(stderr,"istranscript[%d] with abund=%f and path[%d]=%d and nodes[0]=%d and nodes[last]=%d\n",t,transfrag[t]->abundance,i,path[i],transfrag[t]->nodes[0],transfrag[t]->nodes.Last());

					if(i==1 || transfrag[t]->nodes[0]==path[i]) { // first time I encounter transfrag I have to set what abundance to use
						if(!transfrag[t]->real) { // if I still didn't solve transfrag
							transfrag[t]->usepath=-1;
							for(int p=0;p<transfrag[t]->path.Count();p++) {
								int *pos=gpos[edge(transfrag[t]->path[p].node,transfrag[t]->path[p].contnode,gno)];
								if(pos && pathpat[*pos]) {
									transfrag[t]->usepath=p; // this is path dependent
									break;
								}
							}
						}
					}


					if(transfrag[t]->nodes[0]<path[i]) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->abundance;
						if(transfrag[t]->real) capacityleft[i]+=transfrag[t]->abundance;
						else if(transfrag[t]->usepath>-1 && int(transfrag[t]->usepath)<transfrag[t]->path.Count()) { //TODO: this crashes with intv.gtf guides -> to fix
							capacityleft[i]+=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;
						}

						//fprintf(stderr,"add transfrag t=%d i=%d sumleft=%f capacityleft=%f\n",t,i,sumleft[i],capacityleft[i]);
					}
					if(transfrag[t]->nodes.Last()>path[i]) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->abundance;
						if(transfrag[t]->real) capacityright[i]+=transfrag[t]->abundance;
						else if(transfrag[t]->usepath>-1 && int(transfrag[t]->usepath)<transfrag[t]->path.Count())
							capacityright[i]+=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;
						//fprintf(stderr,"add transfrag t=%d i=%d sumright=%f capacityright=%f\n",t,i,sumright[i],capacityright[i]);
					}
				}
				else { // transfrag not on path
					if(path[i]>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(path[i]<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		*/

		if(!capacityleft[i]) return(0);
		if(!capacityright[i]) return(0);

	}

	//if(!full) return(0);

	/*
	{ // DEBUG ONLY
		for(int i=1;i<n-1;i++) {
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	*/

	// compute flow
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		//fprintf(stderr,"i=%d sumright=%f prevflow=%f\n",i,sumright[i],prevflow);
		nodeflux[i]=prevflow/sumright[i];
		//fprintf(stderr,"nodeflux=%f\n",nodeflux[i]);
		capacityright[i]=prevflow;
		prevflow=nodeflux[i]*sumleft[i];
		//fprintf(stderr,"i=%d sumright=%f sumleft=%f prevflow=%f capacityright=%f nodeflux=%f\n",i,sumright[i],sumleft[i],prevflow,capacityright[i],nodeflux[i]);
		//capacityleft[i]=prevflow; // I don't use this
	}

	// * here I don't care what node I treat first
	for(int i=1;i<n-1;i++) if(capacityright[i]){
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {

				float trabundance=transfrag[t]->abundance;
				if(!transfrag[t]->real) {
					if(transfrag[t]->usepath>-1 && int(transfrag[t]->usepath)<transfrag[t]->path.Count())
						trabundance=transfrag[t]->abundance*transfrag[t]->path[int(transfrag[t]->usepath)].abundance;
					else trabundance=0;
				}

				if(trabundance && transfrag[t]->nodes[0]==path[i]) { // transfrag starts at this node
					if(capacityright[i]>trabundance) {
						//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to 0\n",t,transfrag[t]->abundance);
						capacityright[i]-=trabundance;
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=trabundance;
						}
						transfrag[t]->abundance-=trabundance;
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						else if(!transfrag[t]->real) {
							transfrag[t]->path[int(transfrag[t]->usepath)].abundance=0;
							if(transfrag[t]->path.Count()-1 < 2) transfrag[t]->real=true;
							else {
								int np=0;
								for(int p=0;p<transfrag[t]->path.Count();p++)
									if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance) np++;
								if(np<2) transfrag[t]->real=true;
							}
						}
					}
					else {
						//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-capacityright[i]);
						transfrag[t]->abundance-=capacityright[i];
						if(transfrag[t]->abundance<epsilon) {
							transfrag[t]->abundance=0;
						}
						else if(!transfrag[t]->real) {
							//transfrag[t]->path[int(transfrag[t]->usepath)].abundance-=capacityright[i]; // not needed anymore because this stores proportions not actual abundances
							//if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance<epsilon) {
							if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance*transfrag[t]->abundance-capacityright[i]<epsilon) {
								transfrag[t]->path[int(transfrag[t]->usepath)].abundance=0;
								if(transfrag[t]->path.Count()-1 < 2) transfrag[t]->real=true;
								else {
									int np=0;
									for(int p=0;p<transfrag[t]->path.Count();p++)
										if(transfrag[t]->path[int(transfrag[t]->usepath)].abundance) np++;
									if(np<2) transfrag[t]->real=true;
								}
							}
						}

						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=capacityright[i];
						}
						capacityright[i]=0;
						break;
					}
				}


			}
		}
	}

	// I only have to deal with source transfrag
	int nt=no2gnode[path[0]]->trf.Count();
	for(int j=0;j<nt;j++) {
		int t=no2gnode[path[0]]->trf[j];
		if(istranscript[t] && transfrag[t]->abundance) {
			//fprintf(stderr,"Update capacity of transfrag[%d] with value=%f to %f\n",t,transfrag[t]->abundance,transfrag[t]->abundance-prevflow);
			transfrag[t]->abundance-=prevflow;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
			break; // there is no point in updating more than one transfrag from source
		}
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Flow:");
		for(int i=0;i<n;i++)
			fprintf(stderr,"Used %f of node %d[%d]\n",nodeflux[i],i,path[i]);
		fprintf(stderr,"\nTranscript abundances");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	*/

	return(nodeflux[1]);

}



float store_transcript_APPLY_UNISPG(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnodeUnispg>& no2gnode,int& geneno,bool& first,int strand,int gno,GIntHash<int>& gpos, bool& included,
		GBitVec& prevpath, bool full,BundleData *bdata, //float fragno, char* id=NULL) {
		   GffObj* t) {

	fprintf(stderr, ">> Inside store_transcript_APPLY_UNISPG\n");
	float cov=0;
	int len=0;
	CGraphnodeUnispg *prevnode=NULL;
	GVec<GSeg> exons;
	GVec<float> exoncov;
	float excov=0;

	uint refstart=0;
	if(bdata) refstart=(uint)bdata->start;

	// /*
	fprintf(stderr,"store transcript path[0]=%d",path[0]);
	if(t) fprintf(stderr," with id=%s",t->getID());
	fprintf(stderr,"\n");
	// */

	int s=0;
	if(!path[0]) s=1;

	bool firstex=true;
	bool lastex=false;

	for(int i=s;i<path.Count()-1;i++) {
		if(!prevpath[path[i]]) { // if I can find a node that was not included previously in any path then this is a new path
			included=false;
			prevpath[path[i]]=1;
		}
		int *pos=gpos[edge(path[i-1],path[i],gno)]; // if I can find an edge that was not included in any previous path then this is a new path
		if(i && pos && !prevpath[*pos]) {
			included=false;
			if(pos) prevpath[*pos]=1;
		}

		CGraphnodeUnispg *node=no2gnode[path[i]];

	    // moved this one before the usedcov computation since that one wasn't used
		if(t && (node->end<t->start || lastex)) { // I am skipping the nodes that do not overlap the guide so that I don't add them up to coverage
			prevnode=node; continue;
		}

		// push
		float usedcov=nodecov[path[i]]*nodeflux[i]*(node->end-node->start+1);
		//fprintf(stderr,"usedcov=%f for nodecov[path[%d]]=%f nodeflux[%d]=%f node->end=%d node->start=%d\n",usedcov,i,nodecov[path[i]],i,nodeflux[i],node->end,node->start);

		uint nodestart=node->start;
		uint nodeend=node->end;

		if(t) { // I am adjusting the start/end of the exon but shouldn't I also adjust the coverage? -> I added two ifs below
			float firstprop=0;
			float lastprop=0;
			if(firstex) {
				if(nodestart<t->start && bdata) {

					float rightcov=0;

					// cummulative bpcov
					float leftcov=get_cov(2*strand,node->start-refstart,t->start-1-refstart,bdata->bpcov);
					if(node->end>t->start) rightcov=get_cov(2*strand,t->start-refstart,node->end-refstart,bdata->bpcov);

					if(leftcov) firstprop=leftcov/(leftcov+rightcov);
				}
				nodestart=t->start;
			}
			if(t->end<=nodeend || i==path.Count()-2) {
				lastex=true;
				if(t->end<node->end && bdata) {
					//usedcov*=(t->end-nodestart+1)/(nodeend-nodestart+1); // this way I am keeping coverage proportions right
					float leftcov=0;

					// cummulative bpcov
					if(node->start<t->end) {
						leftcov=get_cov(2*strand,node->start-refstart,t->end-refstart,bdata->bpcov);
					}
					float rightcov=get_cov(2*strand,t->end+1-refstart,node->end-refstart,bdata->bpcov);

					if(rightcov) lastprop=rightcov/(leftcov+rightcov);
				}
				nodeend=t->end;
			}
			if(firstprop || lastprop) {
				usedcov-=(firstprop +lastprop)*usedcov;
				nodeflux[i]*=(1-firstprop-lastprop);
			}
		}

		nodecov[path[i]]*=(1-nodeflux[i]); // don't allow this to be less than 0

		if(!prevnode || firstex || node->start>prevnode->end+1) { // this is a new exon
			if(prevnode && !firstex) { // compute exon coverage
				excov/=exons.Last().end-exons.Last().start+1;
				exoncov.Add(excov);
				excov=0;
			}
			GSeg exon(nodestart,nodeend);
			exons.Add(exon);
			firstex=false;
		}
		else if(!firstex) exons.Last().end=nodeend;

		len+=nodeend-nodestart+1;

		cov+=usedcov;
		excov+=usedcov;

		//if(node->cov) fragno+=node->frag*usedcov/node->cov;

		prevnode=node;
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Predicted transcript cov=%f usedcov=%f len=%d path.count=%d ",cov/len, cov,len,path.Count());
		fprintf(stderr,"and exons cov:");
		for(int e=0;e<exons.Count();e++) fprintf(stderr," %d-%d",exons[e].start,exons[e].end);
		fprintf(stderr,"\n");
		if(t) fprintf(stderr,"Ref_id=%s\n",t->getID());
	}
	*/

	// Add last exon coverage
	if(prevnode) { // compute exon coverage
		excov/=exons.Last().end-exons.Last().start+1;
		exoncov.Add(excov);
	}
	if(len) cov/=len;

	//if(t || (cov>=readthr && len>=mintranscriptlen)) { // store transcript here; also accept some coverage fuzziness that would get eliminated later
	//if(t || (cov>=1 && len>=mintranscriptlen)) { // store transcript here; also accept some coverage fuzziness that would get eliminated later
	// sensitive mode:
	if(t || (cov && len>=mintranscriptlen)) { // store transcript here;
		char sign='-';
		if(strand) { sign='+';}
		if(first) { geneno++;}
		//CPrediction *p=new CPrediction(geneno-1, id, exons[0].start, exons.Last().end, cov, sign, fragno, len);
		//if(t) fprintf(stderr,"store prediction with start=%d and end=%d\n",exons[0].start, exons.Last().end);
		float gcov=cov;

		if (t && t->uptr) {
			RC_TData &td = *(RC_TData*) (t->uptr);
			td.in_bundle=3;
			//fprintf(stderr,"st guide %s is stored\n",t->getID());
		}

		/*
		// this was up to version 1.2.1 -> I am not sure about keeping it
		if(t && t->exons.Count()==1) { // if single exon
			RC_TData* tdata=(RC_TData*)(t->uptr);
			if(len) gcov=(tdata->t_exons[0])->movlcount/len;
			if(cov<gcov) gcov=cov;
		}
		*/

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Store transcript in prediction %d with coverage %f \n",pred.Count(),gcov);
			fprintf(stderr,"And exons cov:");
			for(int e=0;e<exons.Count();e++) fprintf(stderr," %g",exoncov[e]);
			fprintf(stderr,"\n");
		}
		*/

		CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, gcov, sign, len);
		p->exons=exons;
		if(t && t->exons.Count()==1) exoncov[0]=gcov;
		p->exoncov=exoncov;
		if(full) p->mergename+='.';
		pred.Add(p);
		first=false;

		//fprintf(stderr,"Transcript stored\n");
	}

	return(cov);
}



// pred[np] = prediction to update; path = the nodes in the prediction or the nodes used from the prediction?
// nodeflux = in store transcript: quantity of transfrags that's used going out from each node;
// nodeflux = here: proportion of the node that is used
// nodecov = coverage of nodes
void update_guide_pred_APPLY_UNISPG(GList<CPrediction>& pred,int np, GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnodeUnispg>& no2gnode,int gno,bool update) {

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Update guide ");
		if(pred[np]->t_eq) fprintf(stderr,"%s ",pred[np]->t_eq->getID());
		fprintf(stderr,"pred[%d]:",np);
		for(int j=0;j<pred[np]->exons.Count();j++) fprintf(stderr," %d-%d",pred[np]->exons[j].start,pred[np]->exons[j].end);
		fprintf(stderr,"\n");
	}
	*/

	int e=0;
	int nex=pred[np]->exons.Count();

	for(int i=0;i<path.Count();i++) {
		if(path[i] && path[i]<gno-1) { // not first or last node in graph
			CGraphnodeUnispg *node=no2gnode[path[i]];
			float addcov=nodeflux[i]*nodecov[path[i]];
			if(update) nodecov[path[i]]-=addcov;
			// addcov/=node->len(); this is not necessary because it was already divided to node length
			while(e<nex && pred[np]->exons[e].end<node->start) e++; // if exon end before start of node skip it
			while(e<nex && pred[np]->exons[e].start<=node->end) { // there is overlap between exon and node
				int ovplen=node->overlapLen(pred[np]->exons[e].start, pred[np]->exons[e].end);
				float excov=addcov*ovplen;
				pred[np]->exoncov[e]+=excov/pred[np]->exons[e].len();
				pred[np]->cov+=excov/abs(pred[np]->tlen);
				e++;
			}
		}
	}
}


float push_guide_maxflow_APPLY_UNISPG(int gno,GVec<int>& path,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,GBitVec& pathpat) {

	float guideabundance=0;

	int n=path.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);
	for(int i=0;i<n;i++) {
		node2path[path[i]]=i;
	}
	GVec<float> capacityleft;   // how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n);
	capacityright.Resize(n);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n);
	sumright.Resize(n);

	/*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"Start push guide max flow algorithm for path ");
		//printBitVec(pathpat);
		fprintf(stderr," :");
		//for(int i=0;i<n;i++) fprintf(stderr," %d:%d",i,path[i]);
		for(int i=0;i<n;i++) fprintf(stderr," %d",path[i]);
		fprintf(stderr,"\n");
		//fprintf(stderr,"Used transcripts:");
		//for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		//fprintf(stderr,"\n");
	}
	*/

	bool marginal=false;

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int nt=no2gnode[path[i]]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[path[i]]->trf[j];
			if(transfrag[t]->abundance) {
				//fprintf(stderr,"Consider transcript %d with abundance %f for node %d\n",t,transfrag[t]->abundance,path[i]);
				if(istranscript[t] || ((pathpat & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;
					//fprintf(stderr,"...on path\n");

					if(marginal && transfrag[t]->nodes[0] && transfrag[t]->nodes.Last()!=gno-1) marginal=false;

					if(transfrag[t]->nodes[0]<path[i]) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->abundance;
						capacityleft[i]+=transfrag[t]->abundance;
					}
					if(transfrag[t]->nodes.Last()>path[i]) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->abundance;
						capacityright[i]+=transfrag[t]->abundance;
					}
				}
				else { // transfrag not on path
					if(path[i]>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(path[i]<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		*/

		if(!capacityleft[i]) return(0);
		if(!capacityright[i]) return(0);


	}

	if(marginal) return(0);

	/*
	{ // DEBUG ONLY
		for(int i=1;i<n-1;i++) {
			fprintf(stderr,"Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
			if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
			else fprintf(stderr,"perc=n/a ");
			fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
			if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
			else fprintf(stderr,"perc=n/a\n");
		}
		fprintf(stderr,"Used transcripts:");
		for(int i=0;i<transfrag.Count();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,transfrag[i]->abundance);
		fprintf(stderr,"\n");
	}
	*/

	// compute flow
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		guideabundance+=no2gnode[path[i]]->cov_s.Get(0)*prevflow/sumright[i];
		prevflow=prevflow*sumleft[i]/sumright[i];
	}

	return(guideabundance);

}



float guidepushflow_APPLY_UNISPG(int g,GVec<CGuide>& guidetrf,int gno,GBitVec& istranscript,GPVec<CTransfrag>& transfrag,
		GPVec<CGraphnodeUnispg>& no2gnode,GVec<float>& nodeflux) {

	int n=guidetrf[g].trf->nodes.Count();
	GVec<int> node2path;
	node2path.Resize(gno,-1);

	//fprintf(stderr,"Process guide %d with pattern: ",g);
	//printBitVec(guidetrf[g].trf->pattern);

	for(int i=0;i<n;i++) {
		node2path[guidetrf[g].trf->nodes[i]]=i;
		nodeflux.cAdd(0.0);
	}

	GVec<float> capacityleft;	// how many transcripts compatible to path enter node
	GVec<float> capacityright;  // how many transcripts compatible to path exit node
	capacityleft.Resize(n);
	capacityright.Resize(n);
	GVec<float> sumleft;        // how many transcripts enter node
	GVec<float> sumright;       // how many transcripts exit node
	sumleft.Resize(n);
	sumright.Resize(n);

	// compute capacities and sums for all nodes
	for(int i=1;i<n-1;i++) {
		int pathi=guidetrf[g].trf->nodes[i];

		int nt=no2gnode[pathi]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[pathi]->trf[j];
			if(transfrag[t]->abundance) {
				if(istranscript[t] || ((guidetrf[g].trf->pattern & transfrag[t]->pattern)==transfrag[t]->pattern)) { // transcript on path
					istranscript[t]=1;
					// check if there are other guides sharing this transcript so that I can allocate proportionally to guide abundances
					float totalcov=guidetrf[g].trf->abundance;
					for(int r=g-1;r>=0;r--) if((guidetrf[r].trf->pattern & transfrag[t]->pattern)==transfrag[t]->pattern) {
						totalcov+=guidetrf[r].trf->abundance;
					}
					float prop=1;
					if(totalcov>guidetrf[g].trf->abundance) prop=guidetrf[g].trf->abundance/totalcov;

					transfrag[t]->usepath=prop*transfrag[t]->abundance;

					if(transfrag[t]->nodes[0]<pathi) { // transfrag starts before this node
						sumleft[i]+=transfrag[t]->usepath;
						capacityleft[i]+=transfrag[t]->usepath;
					}
					if(transfrag[t]->nodes.Last()>pathi) { // transfrag ends after this node
						sumright[i]+=transfrag[t]->usepath;
						capacityright[i]+=transfrag[t]->usepath;
					}
				}
				else { // transfrag not on path
					if(pathi>transfrag[t]->nodes[0]) sumleft[i]+=transfrag[t]->abundance;
					if(pathi<transfrag[t]->nodes.Last()) sumright[i]+=transfrag[t]->abundance;
				}
			}
		}

		if(!capacityleft[i]) return(0);
		if(!capacityright[i]) return(0);
	}

	// compute flow -> deal with rounding errors here
	float prevflow=capacityleft[1];
	for(int i=1;i<n-1;i++) {
		float percleft=prevflow/sumleft[i];
		float percright=capacityright[i]/sumright[i];
		if(percright>percleft) { // more transfrags leave node
			percright=percleft;
		}
		prevflow=percright*sumright[i];
	}
	if(!prevflow) return(0);

	for(int i=n-2;i>0;i--) {
		nodeflux[i]=prevflow/sumright[i];
		if(nodeflux[i]>1) nodeflux[i]=1; // because of rounding errors this could become more than 1
		capacityright[i]=prevflow;
		prevflow=prevflow*sumleft[i]/sumright[i];
		//capacityleft[i]=prevflow;
	}

	for(int i=1;i<n-1;i++) {
		int pathi=guidetrf[g].trf->nodes[i];
		int nt=no2gnode[pathi]->trf.Count();
		for(int j=0;j<nt;j++) {
			int t=no2gnode[pathi]->trf[j];
			if(istranscript[t] && transfrag[t]->abundance) {
				float trabundance=transfrag[t]->usepath;

				if(transfrag[t]->nodes[0]==pathi) { // transfrag starts at this node

					if(capacityright[i]>trabundance) {
						capacityright[i]-=trabundance;
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=trabundance;
						}
						transfrag[t]->abundance-=trabundance;
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						else transfrag[t]->usepath=0; // I only need to delete this if I still have abundance left in transfrag
					}
					else if(capacityright[i]) {
						transfrag[t]->abundance-=capacityright[i];
						if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
						else {
							transfrag[t]->usepath-=capacityright[i];
							if(transfrag[t]->usepath<epsilon)
								transfrag[t]->usepath=0;
						}
						int n2=node2path[transfrag[t]->nodes.Last()];
						for(int k=i+1;k<n2;k++) {
							capacityright[k]-=capacityright[i];
						}
						capacityright[i]=0;
						// break; I can not break because I might still have transfrag that need the cleanup
					}

				}
			}
		}
	}

	// I only have to deal with source transfrag
	int nt=no2gnode[guidetrf[g].trf->nodes[0]]->trf.Count();
	for(int j=0;j<nt;j++) {
		int t=no2gnode[guidetrf[g].trf->nodes[0]]->trf[j];
		if(istranscript[t] && transfrag[t]->abundance) {
			transfrag[t]->abundance-=prevflow;
			if(transfrag[t]->abundance<epsilon) transfrag[t]->abundance=0;
			break;
		}
	} //*/

	return(nodeflux[1]);
}



CTransfrag *find_guide_partial_pat_APPLY_UNISPG(GffObj *guide,GPVec<CGraphnodeUnispg>& no2gnode,int gno,int edgeno,GIntHash<int> &gpos,GVec<int>& olen,int &olensum) {

	CTransfrag *trguide=NULL;
	GBitVec guidepat(gno+edgeno);
	GVec<int> nodes;

	olensum=0;

	int nex=guide->exons.Count();
	int e=0;
	bool terminal=false;
	int lastnode=0; // lastnode intersected by guide
	for(int i=1;i<gno-1;i++) {
		CGraphnodeUnispg *inode=no2gnode[i];
		int len=inode->overlapLen(guide->exons[e]->start,guide->exons[e]->end);
		if(len) { // if node i intersects exon e of guide
			guidepat[i]=1;
			nodes.Add(i);
			olen.Add(len);
			olensum+=len;
			if(lastnode && inode->parentpat[lastnode]) { // the guide has intersected another node in the graph which was i's parent
				if(guide->exons[e]->start==inode->start) {
					if(terminal) { // lastnode ends the same as guide's exon
						int *pos=gpos[edge(lastnode,i,gno)];
						if(pos) guidepat[*pos]=1;
					}
				}
				else if(guide->exons[e]->start<inode->start) {
					if(no2gnode[lastnode]->end+1==inode->start) { // the previous lastnode comes right before this one
						int *pos=gpos[edge(lastnode,i,gno)];
						if(pos) guidepat[*pos]=1;
					}
				}
			}
			lastnode=i;
			terminal=false;
			if(guide->exons[e]->end<=inode->end) {
				if(guide->exons[e]->end==inode->end) terminal=true;
				e++;
				if(e==nex) break; // I've seen all guide's exons
			}
		}
		else if(inode->start>guide->exons[e]->end){ // stay with the same node if node comes after exon
			e++;
			i--;
			if(e==nex) break; // I've seen all guide's exons
		}
	}

	if(olensum) { // guide intersects at least one node in graph
		trguide=new CTransfrag(nodes,guidepat,0,false);
	}

	return(trguide);
}




int store_guide_transcript_APPLY_UNISPG(GList<CPrediction>& pred,GVec<int>& path,GVec<float>& nodeflux,GVec<float>& nodecov,
		GPVec<CGraphnodeUnispg>& no2gnode,int& geneno,bool& first,int gno, GffObj* t,bool update) {

	fprintf(stderr, ">> store_guide_transcript_APPLY_UNISPG! \n");
	// first create the prediction based on the GffObj and then update it's coverage
	GVec<GSeg> exons;
	GVec<float> exoncov;
	int len=0;

	for(int i=0;i<t->exons.Count();i++) {
		exons.Add(t->exons[i]);
		exoncov.cAdd(0.0);
		len+=t->exons[i]->len();
	}

	if(first) { geneno++;first=false;}
	int np=pred.Count();
	CPrediction *p=new CPrediction(geneno-1, t, exons[0].start, exons.Last().end, 0, t->strand, len);
	p->exons=exons;
	p->exoncov=exoncov;
	pred.Add(p);

	if (t && t->uptr) {
		RC_TData &td = *(RC_TData*) (t->uptr);
		td.in_bundle=3;
		//fprintf(stderr,"sg guide %s is stored\n",t->getID());
	}

	update_guide_pred_APPLY_UNISPG(pred,np,path,nodeflux,nodecov,no2gnode,gno,update);

	return(np);
}


bool back_to_source_fast_APPLY_UNISPG(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){

	// fprintf(stderr, ">> Inside back_to_source_fast_APPLY_UNISPG\n");
	// find all parents -> if parent is source then go back
	CGraphnodeUnispg *inode=no2gnode[i];

	int nparents=inode->parent.Count(); // number of parents

	// fprintf(stderr, ">> nparents: %d\n", nparents);

	// /*
	fprintf(stderr,"Parents of node %d are:",i);
	for(int p=0;p<nparents;p++) fprintf(stderr," %d",inode->parent[p]);
	fprintf(stderr,"\n");
	// */

	if(!nparents) {
		// fprintf(stderr, ">> 'back_to_source_fast_APPLY_UNISPG' return true\n");
		return true; // node is source
	}
	int maxp=-1;


	/* I can not do this because then I don't check for continuity
	if(nparents==1) {
		maxp=inode->parent[0];
	}
	else {
	*/
	float maxcov=0;
	//int maxparent=inode->parent[0];
	//float maxparentcov=-1;
	bool exclude=false;
	for(int p=0;p<nparents;p++) {
		float parentcov=0;
		CGraphnodeUnispg *pnode=no2gnode[inode->parent[p]];

		/*
			if(nodecov[inode->parent[p]]>maxparentcov) {
				maxparentcov=nodecov[inode->parent[p]];
				maxparent=inode->parent[p];
			}
		*/

		if(inode->parent[p]==i-1 && i>1 && inode->start==pnode->end+1 &&
				nodecov[i-1]/pnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i-1]) { // adjacent to parent
				// pnode->len() && this is redundant because source should always be 0
				//((nodecov[i-1]/pnode->len() <1000 && nodecov[i]*DROP>nodecov[i-1]) || (pnode->len()<longintronanchor && pnode->parent.Count()==1 && !pnode->parent[0])))  { // adjacent to parent
			exclude=true;
		}
		else {
			pathpat[inode->parent[p]]=1;
			int *pos=gpos[edge(inode->parent[p],i,gno)];
			if(pos) pathpat[*pos]=1;
			else GError("4 Found parent-child %d-%d not linked by edge\n",inode->parent[p],i);

			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				// fprintf(stderr, "Processing j: %d\n", j);
				// fprintf(stderr,"\t >> inside j; parentcov=%f\n",parentcov);
				int t=pnode->trf[j];
				// fprintf(stderr,"\t >> inside j; t=%d\n", t);
				// fprintf(stderr, "\t >> transfrag[t]->abundance: %f\n", transfrag[t]->abundance);
				// fprintf(stderr, "\t >> transfrag[t]->nodes[0]: %d\n", transfrag[t]->nodes[0]);
				// fprintf(stderr, "\t >> inode->parent[p]: %d\n", inode->parent[p]);
				// fprintf(stderr, "\t >> transfrag[t]->nodes.Last(): %d\n", transfrag[t]->nodes.Last());
				// onpath_APPLY_UNISPG(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,inode->parent[p],path[0],no2gnode,gno,gpos);
				// fprintf(stderr, "\t >> After\n");

				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					// fprintf(stderr,"\t >> deleting j:%d\n", j);
					pnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->nodes[0]<=inode->parent[p] && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						//(transfrag[t]->pattern[inode->parent[p]]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath_APPLY_UNISPG(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,inode->parent[p],path[0],no2gnode,gno,gpos)) { // transfrag is compatible with path
					// fprintf(stderr,"\t >> adding abundance at j:%d\n", j);
					// comment: I need to check for the parent to make sure I am not going through another node!!
					//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",inode->parent[p],t,transfrag[t]->abundance);
					parentcov+=transfrag[t]->abundance;
				}
			}

			// fprintf(stderr,"\t >> parentcov: %f;  maxcov: %f\n", parentcov, maxcov);
			if(parentcov>maxcov) {
				maxcov=parentcov;
				maxp=inode->parent[p];
			}
			else if(parentcov==maxcov && maxp!=-1) {
				if(nodecov[maxp]<nodecov[inode->parent[p]]) {
					maxp=inode->parent[p];
				}
			}

			pathpat[inode->parent[p]]=0;
			if(pos) pathpat[*pos]=0;
		}
	}
	if(maxp==-1) {
		//bool goback=false;
		if(exclude && nodecov[i-1]) {
			CGraphnodeUnispg *pnode=no2gnode[i-1];
			float parentcov=0;
			pathpat[i-1]=1;
			int *pos=gpos[edge(i-1,i,gno)];
			if(pos) pathpat[*pos]=1;
			else GError("5 Found parent-child %d-%d not linked by edge\n",i-1,i);

			for(int j=0;j<pnode->trf.Count();j++) { // for all transfrags going through parent
				int t=pnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					pnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->nodes[0]<=i-1 && transfrag[t]->nodes.Last()>=i &&   // transfrag goes from p to i
						//(transfrag[t]->pattern[i-1]||transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath_APPLY_UNISPG(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,i-1,path[0],no2gnode,gno,gpos)) { // transfrag is compatible with path
					// comment: I need to check for the parent to make sure I am not going through another node!!
					//fprintf(stderr,"parent=%d add transfr[%d]->abund=%f\n",i-1,t,transfrag[t]->abundance);
					parentcov+=transfrag[t]->abundance;
				}
			}
			pathpat[i-1]=0;
			if(pos) pathpat[*pos]=0;

			if(parentcov) {
				maxp=i-1;
			}
			else {
				//goback=true;
				// fprintf(stderr,"parentcov: %f\n", parentcov);
				// fprintf(stderr,"return false\n");
				return false;
			}
		}
		else {
			//goback=true;
			// fprintf(stderr,"return false maxp=maxparent=%d\n",maxp);
			return false; //maxp=maxparent;
		}

	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"maxp=%d maxcov=%f\n",maxp,maxcov);
	}
	*/

	if(maxp) { // add maxp to path only if not source
		path.Add(maxp);
	}

	pathpat[maxp]=1;                 // if maxp is source I added it in the pathpat
	int *pos=gpos[edge(maxp,i,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("6 Found parent-child %d-%d not linked by edge\n",maxp,i);

	return back_to_source_fast_APPLY_UNISPG(maxp,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}


bool fwd_to_sink_fast_APPLY_UNISPG(int i,GVec<int>& path,GBitVec& pathpat,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,
		GVec<float>& nodecov,int gno,GIntHash<int>& gpos){
	// fprintf(stderr, ">> Inside fwd_to_sink_fast_APPLY_UNISPG\n");
	// find all parents -> if parent is source then go back
	CGraphnodeUnispg *inode=no2gnode[i];

	int nchildren=inode->child.Count(); // number of children

	/*
	fprintf(stderr,"Children of node %d are:",i);
	for(int c=0;c<nchildren;c++) fprintf(stderr," %d",inode->child[c]);
	fprintf(stderr,"\n");
	*/

	if(!nchildren) return true; // node is sink
	int maxc=-1;

	/* I can not do this because then I don't check for continuity
	if(nchildren==1) {
		maxc=inode->child[0]; // only one child to consider
	}
	else {
	*/
	float maxcov=0;
	//int maxchild=inode->child[0];
	//float maxchildcov=-1;
	bool exclude=false;
	for(int c=0;c<nchildren;c++) {
		float childcov=0;
		CGraphnodeUnispg *cnode=no2gnode[inode->child[c]];
		/*
			if(nodecov[inode->child[c]]>maxchildcov) {
				maxchildcov=nodecov[inode->child[c]];
				maxchild=inode->child[c];
			}
		*/

		if(inode->child[c]==i+1 && i<gno-2 && inode->end+1==cnode->start &&
				nodecov[i+1]/cnode->len() <1000 && nodecov[i]*(DROP+ERROR_PERC)>nodecov[i+1])  { // adjacent to child
				// cnode->len() this is redundant because i<gno-2
				//((nodecov[i+1]/cnode->len() <1000 && nodecov[i]*DROP>nodecov[i+1]) || (cnode->len()<longintronanchor && cnode->child.Count()==1 && cnode->child[0]==gno-1 ) ))  { // adjacent to child
			exclude=true;
		}
		else {
			pathpat[inode->child[c]]=1;
	    	int *pos=gpos[edge(i,inode->child[c],gno)];
	    	if(pos) pathpat[*pos]=1;
	    	else GError("1 Found parent-child: %d-%d not linked by edge!\n",i,inode->child[c]);
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					cnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=inode->child[c] &&   // transfrag goes from i to c
						//(transfrag[t]->pattern[inode->child[c]] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath_APPLY_UNISPG(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],inode->child[c],no2gnode,gno,gpos)) { // transfrag is compatible with path
					//fprintf(stderr,"add transfr[%d]->abund=%f\n",t,transfrag[t]->abundance);
					childcov+=transfrag[t]->abundance;
				}
			}

			if(childcov>maxcov) {
				maxcov=childcov;
				maxc=inode->child[c];
			}
			else if(childcov==maxcov && maxc!=-1) {
				if(nodecov[maxc]<nodecov[inode->child[c]]) {
					maxc=inode->child[c];
				}
			}

			pathpat[inode->child[c]]=0;
			if(pos) pathpat[*pos]=0;
		}
	}
	if(maxc==-1) {
		//bool goback=false;
		if(exclude && nodecov[i+1]) {
			CGraphnodeUnispg *cnode=no2gnode[i+1];
			float childcov=0;
			pathpat[i+1]=1;
	    	int *pos=gpos[edge(i,i+1,gno)];
	    	if(pos) pathpat[*pos]=1;
	    	else GError("2 Found parent-child: %d-%d not linked by edge\n",i,i+1);
			for(int j=0;j<cnode->trf.Count();j++) { // for all transfrags going through child
				int t=cnode->trf[j];
				if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
					cnode->trf.Delete(j);
					j--;
				}
				else if(transfrag[t]->nodes[0]<=i && transfrag[t]->nodes.Last()>=i+1 &&   // transfrag goes from i to c
						//(transfrag[t]->pattern[i+1] || transfrag[t]->pattern[i]) &&    // transfrag is not incomplete through these nodes
						onpath_APPLY_UNISPG(transfrag[t]->pattern,transfrag[t]->nodes,pathpat,path[0],i+1,no2gnode,gno,gpos)) { // transfrag is compatible with path
					childcov+=transfrag[t]->abundance;
				}
			}
			pathpat[i+1]=0; // these were 1 in the previous code, but it doesn't make sense -> I hope it was wrong
			if(pos) pathpat[*pos]=0;

			if(childcov) {
				maxc=i+1;
			}
			else {
				//goback=true;
				//fprintf(stderr,"fwd: return false\n");
				return false;
			}
		}
		else {
			//goback=true;
			//fprintf(stderr,"fwd: return false maxc=maxchild=%d\n",maxc);
			return false; //maxc=maxchild;
		}

		/*if(goback) { // no continuation to source was found
			// check if I have a better path before returning false; might want to restrict to long reads?
			while(path.Count()>1) {
				int c=path.Pop();
				pathpat[c]=0;
				int p=path.Last();
				int *pos=gpos[edge(p,c,gno)];
				if(pos) pathpat[*pos]=0;
				if(no2gnode[p]->child.Last()==gno-1) { // has sink child --> check if there is any abundance left
					for(int j=0;j<no2gnode[p]->trf.Count();j++) { // for all transfrags going through parent
						int t=no2gnode[p]->trf[j];
						if(transfrag[t]->abundance<epsilon) { // this transfrag was used before -> needs to be deleted
							no2gnode[p]->trf.Delete(j);
							j--;
						}
						else if(transfrag[t]->nodes.Last()==gno-1) { // transfrag goes through sink
							pathpat[gno-1]=1;
							int *pos=gpos[edge(p,gno-1,gno)];
							if(pos) pathpat[*pos]=1;
							c=gno-1;
							path.Add(c);
							return true;
						}
					}
				}
			}
			return false;
		}*/

	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"maxc=%d maxcov=%f\n",maxc,maxcov);
	}
	*/

	// add maxp to path
	path.Add(maxc);
	pathpat[maxc]=1;

	int *pos=gpos[edge(i,maxc,gno)];
	if(pos) pathpat[*pos]=1;
	else GError("3 Found parent-child %d-%d not linked by edge\n",i,maxc);

	return fwd_to_sink_fast_APPLY_UNISPG(maxc,path,pathpat,transfrag,no2gnode,nodecov,gno,gpos);
}

bool onpath_APPLY_UNISPG(GBitVec& trpattern,GVec<int>& trnode,GBitVec& pathpattern,int mini,int maxi,GPVec<CGraphnodeUnispg>& no2gnode,int gno,
		GIntHash<int>& gpos) {
	// fprintf(stderr, "Inside onpath_APPLY_UNISPG \n");
	// fprintf(stderr, "trnode[0]: %d\n", trnode[0]);

	if(trnode[0]<mini) // mini can be reached through transcript
	    if(!no2gnode[mini]->parentpat[trnode[0]])	return false;


	if(trnode.Last()>maxi) // from maxi I can reach end of transcript
	    if(!no2gnode[maxi]->childpat[trnode.Last()]) return false;

	int first=1;
	// fprintf(stderr, "After parentpat & childpat\n");

	for(int i=0;i<trnode.Count();i++) {
		if(trnode[i]>=mini && trnode[i]<=maxi) {
	      if(!pathpattern[trnode[i]]) return false;
    	  int *pos=NULL;
    	  if(i) pos=gpos[edge(trnode[i-1],trnode[i],gno)];
	      if(first) {
	    	  first=0;
	    	  if(i && trnode[i]>mini && pos && trpattern[*pos])
	    		  return false;
	    	  if(i && !no2gnode[mini]->parentpat[trnode[i-1]]) return false; // I can not reach mini from previous node in transfrag
	      }
	      else if(i && pos && trpattern[*pos] && !pathpattern[*pos]) return false;
	    }
	    if(trnode[i]>maxi) {
	    	int *pos=gpos[edge(trnode[i-1],trnode[i],gno)];
	    	if(i && trnode[i-1]<maxi && pos && trpattern[*pos])
	    		return false;
	    	if(!no2gnode[maxi]->childpat[trnode[i]]) return false; // I can not reach this node from maxi
	    	break;
	    }
	}

	return true;
}
