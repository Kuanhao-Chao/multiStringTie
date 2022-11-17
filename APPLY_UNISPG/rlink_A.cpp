#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#include "rlink_A.h"


void cov_edge_add_APPLY_UNISPG(GVec<float> *bpcov, int sno, int start, int end, float v) {
	bool neutral=false;
	if(sno!=1) neutral=true; // why is neutral true here: because if the sno is -/+ than I want to add their counts to bpcov[1] too

	// fprintf(stderr, "\t\t>> v: %f\n", v);

	// fprintf(stderr, "\t\tsno: %d; start: %d\n", sno, start);
	// fprintf(stderr, "\t\t	before start => bpcov[sno][start+1]: %f\n", bpcov[sno][start+1]);
	bpcov[sno][start+1]+=v; // if sno==1 then I add v to it here
	// fprintf(stderr, "\t\t	after start => bpcov[sno][start+1]: %f\n", bpcov[sno][start+1]);

	// fprintf(stderr, "\t\tsno: %d; end: %d\n", sno, end);
	// fprintf(stderr, "\t\t	before end => bpcov[sno][end+1]: %f\n", bpcov[sno][end+1]);
	// bpcov[sno][start+1]+=v; // if sno==1 then I add v to it here
	bpcov[sno][end+1]-=v;
	// fprintf(stderr, "\t\t	after end => bpcov[sno][end+1]: %f\n", bpcov[sno][end+1]);


	if(neutral) { // if neutral (i.e. stranded) gets added to bpcov[1] here too => bpcov[1]=all coverage
		bpcov[1][start+1]+=v;
		bpcov[1][end+1]-=v;
	}
}


void add_read_to_cov_APPLY_UNISPG(GList<CReadAln>& rd,int n,GVec<float> *bpcov,int refstart, int refend) {

	// fprintf(stderr, ">>> In add_read_to_cov, (int)rd[n]->strand: %d \n", (int)rd[n]->strand);

	int sno=(int)rd[n]->strand+1; // 0(-),1(.),2(+)

	// fprintf(stderr, ">>> In add_read_to_cov, sno: %d \n", sno);
	float single_count=rd[n]->read_count;
	for(int i=0;i<rd[n]->pair_idx.Count();i++) {
		single_count-=rd[n]->pair_count[i];
		if(rd[n]->pair_idx[i]<=n) { // pair comes before
			int np=rd[n]->pair_idx[i];
			int pcount=rd[np]->segs.Count();
			int rcount=rd[n]->segs.Count();
			int snop=(int)rd[np]->strand+1; //v6
			int p=0;
			int r=0;
			int nsegs=pcount+rcount;
			// fprintf(stderr, ">>> nsegs: %d \n", nsegs);

			int counter = 0;
			while(nsegs) {
				// fprintf(stderr, ">>> loop %d ~~~\n", counter);
				counter++;
				int start;
				int end;
				if(r<rcount) {
					if(p==pcount || rd[np]->segs[p].start>rd[n]->segs[r].end) { // no more pair segments or pair segment comes after read segment
						start=rd[n]->segs[r].start;
						end=rd[n]->segs[r].end;
						r++;
						nsegs--;
					}
					else if(rd[np]->segs[p].end<rd[n]->segs[r].start) {
						start=rd[np]->segs[p].start;
						end=rd[np]->segs[p].end;
						p++;
						nsegs--;
					}
					else { // segments overlap
						start=rd[np]->segs[p].start;
						end=rd[np]->segs[p].end;
						if((int)rd[n]->segs[r].start<start) start=rd[n]->segs[r].start;
						p++;
						nsegs--;
						bool cont=true;
						while(cont) {
							cont=false;
							if(r<rcount && (int)(rd[n]->segs[r].start)<=end) {
								if((int)rd[n]->segs[r].end>end) end=rd[n]->segs[r].end;
								r++;
								nsegs--;
								cont=true;
							}
							if(p<pcount && (int)(rd[np]->segs[p].start)<=end) {
								if((int)rd[np]->segs[p].end>end) end=rd[np]->segs[p].end;
								p++;
								nsegs--;
								cont=true;
							}
						}
					}
				}
				else {
					start=rd[np]->segs[p].start;
					end=rd[np]->segs[p].end;
					p++;
					nsegs--;
				}
				int strand=sno;
				// fprintf(stderr, "\t>>> In add_read_to_cov, sno: %d \n", sno);
				// fprintf(stderr, "\t>>> In add_read_to_cov, snop: %d \n", snop);
				if(sno!=snop) {
					if(sno==1) strand=snop;
					else if(snop!=1) strand=1;
				}
				// fprintf(stderr, "\t>> Iterating pair_idx!!!\n");
				// fprintf(stderr, "\t>>>>> start: %d;  refstart-refend: %d-%d;  start-refstart: %d\n", start, refstart, refend, start-refstart);
				// fprintf(stderr, "\t>>>>> end: %d;  refstart-refend: %d-%d;  end-refstart: %d\n", end, refstart, refend, end-refstart);
				// if (start < refstart) {
				// 	start = refstart;
				// }
				// if (end > refend) {
				// 	end = refend-1;
				// }

        		// if (start >= refstart && end <= refend) {
					cov_edge_add_APPLY_UNISPG(bpcov,strand,start-refstart,end-refstart+1,rd[n]->pair_count[i]);
				// }
			}
		}
	}
	if(single_count>epsilon) {
		// fprintf(stderr, ">> single_count>epsilon\n");
		for(int i=0;i<rd[n]->segs.Count();i++) {

			if (int(rd[n]->segs[i].end) < refstart) continue;				

			int start = int(rd[n]->segs[i].start);
			int end = int(rd[n]->segs[i].end);
			// fprintf(stderr, "\t\tstart: %d;  end: %d\n", start, end);
			// fprintf(stderr, "\t\trefstart: %d;  refend: %d\n", refstart, refend);



			if (start <= refstart && end <= refstart) {
				/*********************************************
				** 1. (r)---------(r)   |(ref).............(ref)|
				**********************************************/
				// Do nothing
				// fprintf(stderr, ">>> (1) Do not process the nodes\n");
			} else if (start >= refend && end >= refend) {
				/*********************************************
				** 1. |(ref).............(ref)|  (r)---------(r)
				**********************************************/
				// Do nothing
				// fprintf(stderr, ">>> (2) Do not process the nodes\n");
				// fprintf(stderr, "@@@@ read: %d - %d\n", rd[n]->start, rd[n]->end);
			} else {
				if (start < refstart) {
					// fprintf(stderr, "!!!! Warning!!!\n");
					start = refstart;
				}
				if (end > refend) {
					end = refend-1;
				}
				// fprintf(stderr, ">>>>> start: %u;  refstart: %d;  start-refstart: %d; refend: %d\n", rd[n]->segs[i].start, refstart, rd[n]->segs[i].start-refstart, refend);
				// fprintf(stderr, ">>>>> end: %u;  refstart: %d;  end-refstart: %d; refend: %d\n", rd[n]->segs[i].end, refstart, rd[n]->segs[i].end-refstart, refend);
				cov_edge_add_APPLY_UNISPG(bpcov, sno, start-refstart, end-refstart+1, single_count);
			}
		}
	}
}




// compute coverages and junction support here
void count_good_junctions_APPLY_UNISPG(BundleData* bdata) {

	GList<CReadAln>& readlist = bdata->readlist;
	GList<CJunction>& junction = bdata->junction;
	GVec<float>* bpcov = bdata->bpcov;
	int refstart=bdata->start;
	int refend=bdata->end;
	bool modified=false;
	//GHash<CJunction*> jhash(false);
	GHashMap<CJunction*, CJunction*> jhash(false); //hash of pointers
	//char sbuf[20];


	for(int s=0;s<3;s++) bpcov[s].Resize(refend-refstart+3);

	GVec<int> unstranded; // remembers unstranded reads

	/*****************************
	 * Iterating all reads
	 *****************************/
	for(int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		// fprintf(stderr, ">> rd[n]->pair_idx.Count(): %d\n", rd.pair_idx.Count());
		// fprintf(stderr, ">> n: %d\n", n);

		float rdcount=rd.read_count;

		int nex=rd.segs.Count();

		/*****************************
		 * Process the read only if it is "completely" inside the bundle.
		 *****************************/
        // if (rd.start >= refstart && rd.end <= refend) {
			add_read_to_cov_APPLY_UNISPG(readlist,n,bpcov,refstart, refend);
			if(rdcount>1) rdcount=1;

			GVec<uint> leftsup;
			GVec<uint> rightsup;
			uint maxleftsupport=0;
			uint maxrightsupport=0;

			//int sno=(int)rd.strand+1; // 0(-),1(.),2(+)
			//if(nex>1) fprintf(stderr,"Process spliced read[%d] with cov=%f and sno=%d: ",n,rdcount,sno);
			for(int i=0;i<nex;i++) {
				//if(nex>1) fprintf(stderr," %d-%d",rd.segs[i].start,rd.segs[i].end);
				if(i) {
					//fprintf(stderr,":%d",rd.juncs[i-1]->strand);

					if(modified) { // see if read uses modified junction -> correct it
						//sprintf(sbuf, "%p", rd.juncs[i-1]);
						CJunction* jp=jhash[rd.juncs[i-1]];
						if(jp) {
							if(rd.segs[i-1].start>jp->start || rd.segs[i].end<jp->end) {

								if(rd.segs[i-1].start<=jp->start) rd.segs[i-1].end=jp->start;
								if(rd.segs[i].end>=jp->end) rd.segs[i].start=jp->end;

								rd.juncs[i-1]->nreads-=rd.read_count;
								//if(rd.juncs[i-1]->nreads<ERROR_PERC) rd.juncs[i-1]->strand=0; // this approach removes a perfectly valid (guide) junction which might not be the desired effect
								rd.juncs[i-1]=junction[0];

							}
							else {
								rd.juncs[i-1]=jp;
								if(!rd.strand) rd.strand=jp->strand;
								//fprintf(stderr," [correct rd from %d-%d to %d-%d]",rd.segs[i-1].end,rd.segs[i].start,jp->start,jp->end);
								if(rd.segs[i-1].start<=jp->start) rd.segs[i-1].end=jp->start;
								if(rd.segs[i].end>=jp->end) rd.segs[i].start=jp->end;
							}

						}
						else {
							if(rd.segs[i-1].start>rd.juncs[i-1]->start || rd.segs[i].end<rd.juncs[i-1]->end) {
								if(rd.segs[i-1].end!=rd.juncs[i-1]->start && rd.segs[i-1].start<=rd.juncs[i-1]->start) rd.segs[i-1].end=rd.juncs[i-1]->start;
								if(rd.segs[i].start!=rd.juncs[i-1]->end && rd.segs[i].end>=rd.juncs[i-1]->end) rd.segs[i].start=rd.juncs[i-1]->end;
								rd.juncs[i-1]->nreads-=rd.read_count;
								//if(rd.juncs[i-1]->nreads<ERROR_PERC) rd.juncs[i-1]->strand=0; // this approach removes a perfectly valid (guide) junction which might not be the desired effect
								rd.juncs[i-1]=junction[0];
							}
							else {
								//if(rd.segs[i-1].end!=rd.juncs[i-1]->start || rd.segs[i].start!=rd.juncs[i-1]->end) fprintf(stderr," [chg rd from %d-%d to %d-%d]",rd.segs[i-1].end,rd.segs[i].start,rd.juncs[i-1]->start,rd.juncs[i-1]->end);
								if(rd.segs[i-1].end!=rd.juncs[i-1]->start && rd.segs[i-1].start<=rd.juncs[i-1]->start) rd.segs[i-1].end=rd.juncs[i-1]->start;
								if(rd.segs[i].start!=rd.juncs[i-1]->end && rd.segs[i].end>=rd.juncs[i-1]->end) rd.segs[i].start=rd.juncs[i-1]->end;
							}
						}
					}

					if(rd.segs[i-1].len()>maxleftsupport) maxleftsupport=rd.segs[i-1].len();
					if(rd.segs[nex-i].len()>maxrightsupport) maxrightsupport=rd.segs[nex-i].len();
					leftsup.Add(maxleftsupport);
					rightsup.Add(maxrightsupport);
					//rd.juncs[i-1]->nreads+=rd.read_count;
					//if(rd.unitig) rd.juncs[i-1]->guide_match=true; // v7 this might be a little too much!
				}
				//if(!rd.unitig) cov_edge_add_APPLY_UNISPG(bpcov,sno,rd.segs[i].start-refstart,rd.segs[i].end+1-refstart,rd.read_count);

			}
			//if(nex>1) fprintf(stderr," With anchors: ");
			for(int i=1;i<nex;i++) {
				//if(!rd.juncs[i-1]->strand) { rd.juncs[i-1]->strand = rd.strand; }
				uint anchor=junctionsupport;
				//if(rd.juncs[i-1]->len()>longintron || rd.juncs[i-1]->nreads-rd.juncs[i-1]->mm<junctionthr || rd.juncs[i-1]->nreads-rd.juncs[i-1]->nm<junctionthr) {
				if(rd.juncs[i-1]->len()>longintron || (rd.juncs[i-1]->nreads-rd.juncs[i-1]->nm<junctionthr && !rd.juncs[i-1]->mm)) {
					//if(rd.juncs[i-1]->len()>verylongintron) { if(anchor<verylongintronanchor) anchor=verylongintronanchor; } // I want to use a longer anchor for long introns to believe them
					//else
						if(anchor<longintronanchor) anchor=longintronanchor;
				}

				//if(rd.unitig) anchor=1; // v7 also **** if not trimming involved in super-read creation then comment this; unitigs should be trimmed so I will accept them as anchors

				//if(leftsup[i-1]>=anchor && rightsup[nex-i-1]>=anchor) rd.juncs[i-1]->nreads_good+=rd.read_count;
				if(leftsup[i-1]>=anchor) { // support only comes from spliced reads that are bigger than the anchor
					rd.juncs[i-1]->leftsupport+=rdcount;
					if(rightsup[nex-i-1]>=anchor) {
						rd.juncs[i-1]->rightsupport+=rdcount;
						rd.juncs[i-1]->nreads_good+=rdcount;
						//rd.juncs[i-1]->nreads_good+=rd.read_count*rd.nh; // try to see if this helps
						//if(leftsup[i-1]>=anchor+1/ERROR_PERC && rightsup[nex-i-1]>=anchor+1/ERROR_PERC) rd.juncs[i-1]->strong=true;
					}
				}
				else if(rightsup[nex-i-1]>=anchor) {
					rd.juncs[i-1]->rightsupport+=rdcount;
				}
				//fprintf(stderr," %d(%d-%d)[%f-%f][%d][%f][%f]",anchor,leftsup[i-1],rightsup[nex-i-1],rd.juncs[i-1]->leftsupport,rd.juncs[i-1]->rightsupport,rd.juncs[i-1]->strand,rd.juncs[i-1]->nm,rd.juncs[i-1]->nreads);
			}
			//if(nex>1) fprintf(stderr,"\n");
		// }
	}


	//fprintf(stderr,"Compute coverages:\n");

	// this code enssures that I can find interval coverages very fast
	int m=int((bpcov[1].Count()-1)/BSIZE);
	//fprintf(stderr,"m=%d refstart=%d refend=%d bpcount=%d\n",m,refstart,refend,bpcov[1].Count());
	int k=0;
	float prev_val[3]={0,0,0};
	float prev_sum[3]={0,0,0};
	while(k<m) {
		int end=(k+1)*BSIZE;

		//fprintf(stderr,"end=%d\n",end);
		for(int i=k*BSIZE;i<end;i++) {

			for(int s=0;s<3;s++) {
				//fprintf(stderr,"(1)bpcov[%d][%d]=%f ",s,i,bpcov[s][i]);
				bpcov[s][i]+=prev_val[s];
				if(bpcov[s][i]<0) bpcov[s][i]=0;
				prev_val[s]=bpcov[s][i];
				//fprintf(stderr,"(2)bpcov[%d][%d]=%f ",s,i,bpcov[s][i]);
				bpcov[s][i]+=prev_sum[s];
				prev_sum[s]=bpcov[s][i];
				//fprintf(stderr,"(2)bpcov[%d][%d]=%f ",s,i,bpcov[s][i]);
				//fprintf(stderr,"bpPos[%d][%d]: prev_val=%f prev_sum=%f\n",s,i,prev_val[s],prev_sum[s]);
			}
		}
		for(int s=0;s<3;s++) prev_sum[s]=0;
		k++;
	}

	for(int i=k*BSIZE;i<bpcov[1].Count();i++)
		for(int s=0;s<3;s++) {
			bpcov[s][i]+=prev_val[s];
			if(bpcov[s][i]<0) bpcov[s][i]=0;
			prev_val[s]=bpcov[s][i];
			bpcov[s][i]+=prev_sum[s];
			prev_sum[s]=bpcov[s][i];
			//fprintf(stderr,"bpPos[%d][%d]: cov=%f sumbpcov=%f\n",s,i,prev_val[s],prev_sum[s]);
		}
}


CTransfrag *findtrf_in_treepat_APPLY_UNISPG(int gno,GIntHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) { // doesn't work for patterns including source node
	// fprintf(stderr, ">> Inside 'findtrf_in_treepat_APPLY_UNISPG'\n");
	return NULL;
	if(!tr2no) return(NULL);

	CTreePat *tree=tr2no;
	for(int n=0;n<node.Count();n++) {
		if(n) { // not the first node in pattern
			int *pos=gpos[edge(node[n-1],node[n],gno)];
			if(pos && pattern[*pos]) { // there is an edge between node[n-1] and node[n]
				tree=tree->nextpat[gno-1-node[n-1]+node[n]-node[n-1]-1];
			}
			else tree=tree->nextpat[node[n]-node[n-1]-1];
		}
		else tree=tree->nextpat[node[n]-1];
		if(!tree) return(NULL);
	}
	return(tree->tr);
}