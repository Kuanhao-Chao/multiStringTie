#include "build_graphs_C.h"

/*****************************
 * CREATE_UNISPG
 *****************************/
int build_graphs_CREATE_UNISPG(BundleData* bdata, UnispgGp_CREATE* unispg_gp, int fidx) {
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
	GVec<int> *readgroup = new GVec<int>[readlist.Count()]; // remebers groups for each read; don't forget to delete it when no longer needed
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
			for(int i=0;i<tdata->t_introns.Count();i++) {
				if(!tdata->t_introns[i]->rcount) {
					covered=false;
					break;
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
		// fprintf(stderr,"check junction:%d-%d:%d leftsupport=%f rightsupport=%f nm=%f nreads=%f\n",junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->leftsupport,junction[i]->rightsupport,junction[i]->nm,junction[i]->nreads);

		if((!higherr) && junction[i]->strand && junction[i]->nm==junction[i]->nreads && !junction[i]->guide_match) {
			higherr=true;
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
				// fprintf(stderr,"junction:%d-%d:%d is %c%c-%c%c\n",junction[i]->start,junction[i]->end,junction[i]->strand,bdata->gseq[junction[i]->start+1-refstart],bdata->gseq[junction[i]->start+2-refstart],bdata->gseq[junction[i]->end-2-refstart],bdata->gseq[junction[i]->end-1-refstart]);
			}
		}
		else if(junction[i]->strand) leftsupport[(1+junction[i]->strand)/2]+=junction[i]->leftsupport;
		// fprintf(stderr,"leftsupport[%d]=%f\n",(1+junction[i]->strand)/2,leftsupport[(1+junction[i]->strand)/2]);


		// fprintf(stderr,"check ejunction:%d-%d:%d leftsupport=%f rightsupport=%f nm=%f nreads=%f\n",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->leftsupport,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);
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
				// fprintf(stderr,"ejunction:%d-%d:%d is %c%c-%c%c\n",ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,bdata->gseq[ejunction[i]->start+1-refstart],bdata->gseq[ejunction[i]->start+2-refstart],bdata->gseq[ejunction[i]->end-2-refstart],bdata->gseq[ejunction[i]->end-1-refstart]);
			}

		}
		else if(ejunction[i]->strand) rightsupport[(1+ejunction[i]->strand)/2]+=ejunction[i]->rightsupport;
		// fprintf(stderr,"rightsupport[%d]=%f\n",(1+ejunction[i]->strand)/2,rightsupport[(1+ejunction[i]->strand)/2]);
	}
	// end adjusting leftsupport and rightsupport

	// fprintf(stderr,"junction support computed\n");

	/*****************************
	 ** Step 2: there are some reads that contain very bad junctions -> need to find better closest junctions
	 *****************************/
	if(higherr) { 
		uint juncsupport=junctionsupport;
		//fprintf(stderr,"In higherr!\n");
		GVec<int> jstarts; // good starting junctions
		GVec<int> jends; // good ending junctions

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
			if(changeright) {
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
						if(ejunction[ek]->nreads_good<0) {
							newend=rd.segs[i+1].end+1;
						}
						else {
							newend=ejunction[ek]->end;
							//fprintf(stderr,"junction has newend=%d from junction[ek=%d]\n",newend,ek);
							if(ejunction[ek]->strand) ek=-1;
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
					// fprintf(stderr, "read[%d] adjusted to junction:%d-%d\n",n,rd.segs[i].end,rd.segs[i+1].start);
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


		// if(rd.juncs.Count()) fprintf(stderr,"read[%d] keep=%d\n",n,keep);
		// if(rd.strand) fprintf(stderr,"read[%d] has strand %d\n",n,rd.strand);


		//if(keep) { // if it's a good read that needs to be kept


			/*
			fprintf(stderr,"add read %d:%d-%d w/count=%g for color=%d with npairs=%d\n",n,readlist[n]->start,readlist[n]->end,readlist[n]->read_count,color,readlist[n]->pair_idx.Count());
			fprintf(stderr,"add read[%d]:%d-%d:%d w/count=%g w/exons:",n,readlist[n]->start,readlist[n]->end,readlist[n]->strand,readlist[n]->read_count);
			for(i=0;i<rd.juncs.Count();i++) { fprintf(stderr," %d-%d:%d",rd.segs[i].start,rd.segs[i].end,rd.juncs[i]->strand);}
			fprintf(stderr," %d-%d\n",rd.segs[i].start,rd.segs[i].end);
			*/

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

			// fprintf(stderr,"now color=%d\n",color);
		//}
		//else { fprintf(stderr,"read[%d] is not kept\n",n);}
		//else clean_read_junctions(readlist[n]);
	}

	if(resort) {
		junction.setSorted(true);
		ejunction.setSorted(juncCmpEnd);
	}
	// fprintf(stderr,"fragno=%d fraglen=%g\n",fragno,fraglen);
	//if(fragno) fraglen/=fragno;


	/*****************************
	 ** Step 4: 'merge_fwd_groups' function
	 ** 	merge groups that are close together or __groups that are within the same exon of a reference gene__
	 *****************************/
	if(bundledist || (guides.Count())) {
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

		// fprintf(stderr, ">> grcol before: %d\n", grcol);			fprintf(stderr, ">> eqposcol count: %d\n", eqposcol.Count());
		while(equalcolor[grcol]!=grcol) {
			grcol=equalcolor[grcol];
			// fprintf(stderr, ">> grcol in while: %d\n", grcol);
		}
		// fprintf(stderr, ">> grcol after: %d\n", grcol);

		// while(equalcolor[grcol]!=grcol) {
		// 	grcol=equalcolor[grcol];
		// }
		currgroup[nextgr]->color=grcol;

		// fprintf(stderr,"group %d id=%d: %u-%u col=%d from col=%d\n",nextgr,currgroup[nextgr]->grid,currgroup[nextgr]->start,currgroup[nextgr]->end,grcol,prevcol);

		if(nextgr == 1) { // unknown strand group

			if(prevgroup[0]!=NULL && currgroup[nextgr]->start <= prevgroup[0]->end+bundledist) { // overlaps previous negative group ; this needs bundledist
				// fprintf(stderr,"\tovlp to neg group: %u-%u\n",prevgroup[0]->start,prevgroup[0]->end);
				set_strandcol(currgroup[nextgr],prevgroup[0],prevgroup[0]->color,eqnegcol,equalcolor);
				uint maxstart = currgroup[nextgr]->start > prevgroup[0]->start ? currgroup[nextgr]->start : prevgroup[0]->start;
				uint minend = currgroup[nextgr]->end < prevgroup[0]->end ? currgroup[nextgr]->end : prevgroup[0]->end;
				if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0

				// fprintf(stderr, "## prevgroup[0]->cov_sum: %f\n", prevgroup[0]->cov_sum);
				// fprintf(stderr, "## currgroup[nextgr]->neg_prop: %f\n", currgroup[nextgr]->neg_prop);
				// fprintf(stderr, "## prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len(): %f\n", prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len());

				currgroup[nextgr]->neg_prop+=prevgroup[0]->cov_sum*(minend-maxstart+1)/prevgroup[0]->len();
			}

			while(currgroup[0]!=NULL && currgroup[nextgr]->start <= currgroup[0]->end+bundledist && currgroup[0]->start <= currgroup[nextgr]->end +bundledist) { // overlaps current negative strand group
				// fprintf(stderr,"\tovlp to neg group: %u-%u\n",currgroup[0]->start,currgroup[0]->end);

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

				// fprintf(stderr, "## currgroup[0]->cov_sum: %f\n", currgroup[0]->cov_sum);
				// fprintf(stderr, "## currgroup[nextgr]->neg_prop: %f\n", currgroup[nextgr]->neg_prop);
				// fprintf(stderr, "## currgroup[0]->cov_sum*(minend-maxstart+1)/currgroup[0]->len(): %f\n", currgroup[0]->cov_sum*(minend-maxstart+1)/currgroup[0]->len());


				currgroup[nextgr]->neg_prop+=currgroup[0]->cov_sum*(minend-maxstart+1)/currgroup[0]->len();


				prevgroup[0]=currgroup[0];
				currgroup[0]=currgroup[0]->next_gr;
			}

			float pos_prop=0;
			if(prevgroup[2]!=NULL && currgroup[nextgr]->start <= prevgroup[2]->end + bundledist) { // overlaps positive strand group
				// fprintf(stderr,"\tovlp to pos group: %u-%u\n",prevgroup[2]->start,prevgroup[2]->end);
				set_strandcol(currgroup[nextgr],prevgroup[2],prevgroup[2]->color,eqposcol,equalcolor);
				if(currgroup[nextgr]->neg_prop) {
					uint maxstart = currgroup[nextgr]->start > prevgroup[2]->start ? currgroup[nextgr]->start : prevgroup[2]->start;
					uint minend = currgroup[nextgr]->end < prevgroup[2]->end ? currgroup[nextgr]->end : prevgroup[2]->end;
					if(minend<maxstart) minend=maxstart; // this can only happen if bundledist >0
					pos_prop+=prevgroup[2]->cov_sum*(minend-maxstart+1)/prevgroup[2]->len();
				}
			}

			while(currgroup[2]!=NULL && currgroup[nextgr]->start <= currgroup[2]->end +bundledist && currgroup[2]->start <= currgroup[nextgr]->end + bundledist) { // overlaps positive strand group
				// fprintf(stderr,"\tovlp to pos group: %u-%u\n",currgroup[2]->start,currgroup[2]->end);

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
			// fprintf(stderr,"neg_prop=%g pos_prop=%g\n",currgroup[nextgr]->neg_prop,pos_prop);
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

	/*
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
	*/
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

		if(bundle[1][b]->cov && ((bundle[1][b]->multi/bundle[1][b]->cov)<=mcov*(1-ERROR_PERC) || guide_ovlp)) { // && (guides.Count() || adaptive || bundle[1][b]->len >= mintranscriptlen)) { // there might be small transfrags that are worth showing, but here I am ignoring them
    		// bundle might contain multiple fragments of a transcript but since we don't know the complete structure -> print only the pieces that are well represented
    		CBundlenode *currbnode=bnode[1][bundle[1][b]->startnode];
    		int t=1;
			// fprintf(node_unispg_unstrand_bed,"strict digraph %d_%d {", refstart, refend);
    		while(currbnode!=NULL) {

				/*****************************
				 ** 3. Write out global splice graph in DOT format
				*****************************/
				/****************
				 **  KH Adding 
				****************/
				// graphno[s][b]: number of nodes in graph.
                // fprintf(node_unispg_unstrand_bed, "chr22\t%d\t%d\t%d\t%f\t%s\n", currbnode->start,currbnode->end, currbnode->bid, currbnode->cov, ".");

				// fprintf(node_unispg_unstrand_bed,"%d[start=%d end=%d cov=%f];",currbnode->bid,currbnode->start,currbnode->end,currbnode->cov);
				/****************
				 **  END KH Adding 
				****************/

    			//int len=currbnode->end-currbnode->start+1;
    			//float cov=currbnode->cov/(currbnode->end-currbnode->start+1);
    			bool printguides=false;

				for(int i=0;i<bnodeguides[currbnode->bid].Count();i++) {
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
									p->exoncov.Add(cov);
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
									p->exoncov.Add(cov);
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
						p->exoncov.Add(cov);
    					pred.Add(p);
    					t++;
    				}
    			}
    			currbnode=currbnode->nextnode;
    		}
			// fprintf(node_unispg_unstrand_bed,"}\n");
    	}
    }
    //fprintf(stderr,"Done with unstranded bundles\n");
    if (bnodeguides) delete[] bnodeguides;



	/*****************************
	 ** Step 11: Defining parameters here!!!!
	 ** 	build graphs for stranded bundles here
	 *****************************/
    if(startgroup[0]!=NULL || startgroup[2]!=NULL) {  // Condition 1: there are stranded groups to process
		/*****************************
		 ** I don't need the groups here anymore : I only use their numbers
		*****************************/
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
			// fprintf(stderr, "** graph_no: %d\n", graph_no);
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

    				// fprintf(stderr,"Bundle is: %d - %d start at g=%d sno=%d b=%d\n",bnode[sno][bundle[sno][b]->startnode]->start,bnode[sno][bundle[sno][b]->lastnodeid]->end,g,sno,b);
					// fprintf(stderr, "bnode[%d][bundle[%d][%d]->startnode]->start: %d \n", sno, sno, b, bnode[sno][bundle[sno][b]->startnode]->start);
					// fprintf(stderr, "bnode[%d][bundle[%d][%d]->lastnodeid]->end: %d \n", sno, sno, b, bnode[sno][bundle[sno][b]->lastnodeid]->end);

    				while(g<ng && guides[g]->end<bnode[sno][bundle[sno][b]->startnode]->start) g++;

    				int cg=g;
    				int nolap=0;

					/****************
					 **  This is reference guide. Ignore first.
					 ****************/
    				while(cg<ng && guides[cg]->start<=bnode[sno][bundle[sno][b]->lastnodeid]->end) { // this are potential guides that might overlap the current bundle, and they might introduce extra edges
    					// fprintf(stderr,"...consider guide cg=%d with strand=%c and in_bundle=%d\n",cg,guides[cg]->strand,((RC_TData*)(guides[cg]->uptr))->in_bundle);
    					if((guides[cg]->strand==strnd || guides[cg]->strand=='.') && ((RC_TData*)(guides[cg]->uptr))->in_bundle>=2) {
    						// fprintf(stderr,"Add guide g=%d with start=%d end=%d\n",cg,guides[cg]->start,guides[cg]->end);
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

						// no2gnode / transfrag / gpos
    					graphno[s][b]=create_graph_unispg(refstart,s,b,bundle[sno][b],bnode[sno],junction,ejunction,
    							bundle2graph,no2gnode,transfrag,gpos,bdata,edgeno[s][b],lastgpos[s][b],guideedge); // also I need to remember graph coverages somewhere -> probably in the create_graph procedure

						// tr2no
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
    	// fprintf(stderr,"Done creating graphs\n");
    	/*
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
    	*/
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
				// fprintf(stderr, ">> single_count: %f\n", single_count);
    			for(int j=0; j<readlist[n]->pair_idx.Count();j++) {
    				int np=readlist[n]->pair_idx[j];
    				if(np>-1) {
    					single_count-=readlist[n]->pair_count[j];
    					if(n<np) {
							// no2gnode / transfrag / gpos / tr2no
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
		 ** 3-1. Creating the Unispg for the universal graph and add it into UnispgGp
		*****************************/
        if (multiMode) {
            for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> thoses shouldn't have junctions
                int s=sno/2; // adjusted strand due to ignoring neutral strand				
				/*****************************
				 ** 3-2. Write out global splice graph in DOT format
				*****************************/
				unispg_gp->AddGraph(fidx, s, bdata->refseq, no2gnode[s], bundle[sno].Count());

				// unispg_gp->construct_transfrag_unispg(fidx, s);

				// unispg_gp->construct_treepat_unispg(graphno[s][b],gpos[s][b],transfrag[s][b]);
    			// unispg_gp->get_fragment_pattern(readlist,n,-1,single_count,readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);

				// create_graph_unispg
				// 1. construct_treepat_unispg
				// 2. get_fragment_pattern
            }

			unispg_gp->WriteUNISPG_DOT(fidx, bdata->refseq);
			unispg_gp->graph_num[0] = 0;
			unispg_gp->graph_num[1] = 0;
			unispg_gp->Clear_no2gnode_unispg();
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

					// no2gnode / transfrag / tr2no / gpos
    				process_transfrags_unispg(s,graphno[s][b],edgeno[s][b],no2gnode[s][b],transfrag[s][b],tr2no[s][b],gpos[s][b],guidetrf,pred,trflong);
    				//get_trf_long_unispg(graphno[s][b],edgeno[s][b], gpos[s][b],no2gnode[s][b],transfrag[s][b],geneno,s,pred,trflong);


    				/*
    				{ //DEBUG ONLY
    					//printTime(stderr);
    					fprintf(stderr,">>> There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
    					for(int i=1;i<graphno[s][b]-1;i++) {
    						
							fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f capacity=%f rate=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len(), no2gnode[s][b][i]->capacity, no2gnode[s][b][i]->rate);

    						fprintf(stderr," parents:");
    						for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
    						fprintf(stderr," children:");
    						for(int j=0;j<no2gnode[s][b][i]->child.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->child[j]);
    						fprintf(stderr," trf=");
    						for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
    						fprintf(stderr,"\n");

							float node_coverage_neg = get_cov(0, no2gnode[s][b][i]->start-refstart, no2gnode[s][b][i]->end-refstart, bdata->bpcov);
							float node_coverage_uns = get_cov(1, no2gnode[s][b][i]->start-refstart, no2gnode[s][b][i]->end-refstart, bdata->bpcov);
							float node_coverage_pos = get_cov(2, no2gnode[s][b][i]->start-refstart, no2gnode[s][b][i]->end-refstart, bdata->bpcov);
							node_coverage_uns = node_coverage_uns - node_coverage_neg - node_coverage_pos;
							float node_coverage_neg_sign = get_cov(0, no2gnode[s][b][i]->start-refstart, no2gnode[s][b][i]->end-refstart, bdata->bpcov);
							float node_coverage_pos_sign = get_cov_sign(2, no2gnode[s][b][i]->start-refstart, no2gnode[s][b][i]->end-refstart, bdata->bpcov);


    						fprintf(stderr,"\t && node_coverage_neg=%f  node_coverage_uns=%f  node_coverage_pos=%f\n", node_coverage_neg, node_coverage_uns, node_coverage_pos);
    						fprintf(stderr,"\t && node_coverage_neg_sign=%f  node_coverage_pos_sign=%f\n", node_coverage_neg_sign, node_coverage_pos_sign);
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
    				*/

/*
#ifdef GMEMTRACE
    				double vm,rsm;
    				get_mem_usage(vm, rsm);
    				GMessage("\t\tM(after process_transfrags):build_graphs memory usage: rsm=%6.1fMB vm=%6.1fMB\n",rsm/1024,vm/1024);
#endif
*/

    				// fprintf(stderr,"guidetrf no=%d\n",guidetrf.Count());

    				// find transcripts now

					// gpos / no2gnode / transfrag
					geneno=find_transcripts_unispg(graphno[s][b],edgeno[s][b],gpos[s][b],no2gnode[s][b],transfrag[s][b],
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
