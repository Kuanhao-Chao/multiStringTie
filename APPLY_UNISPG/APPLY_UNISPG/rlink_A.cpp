#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#include "rlink_A.h"


void cov_edge_add_APPLY_UNISPG(GVec<float> *bpcov, int sno, int start, int end, float v) {
	bool neutral=false;
	if(sno!=1) neutral=true; // why is neutral true here: because if the sno is -/+ than I want to add their counts to bpcov[1] too

	fprintf(stderr, "\t\t>> v: %f\n", v);

	fprintf(stderr, "\t\tsno: %d; start: %d\n", sno, start);
	fprintf(stderr, "\t\t	before start => bpcov[sno][start+1]: %f\n", bpcov[sno][start+1]);
	bpcov[sno][start+1]+=v; // if sno==1 then I add v to it here
	fprintf(stderr, "\t\t	after start => bpcov[sno][start+1]: %f\n", bpcov[sno][start+1]);

	fprintf(stderr, "\t\tsno: %d; end: %d\n", sno, end);
	fprintf(stderr, "\t\t	before end => bpcov[sno][end+1]: %f\n", bpcov[sno][end+1]);
	// bpcov[sno][start+1]+=v; // if sno==1 then I add v to it here
	bpcov[sno][end+1]-=v;
	fprintf(stderr, "\t\t	after end => bpcov[sno][end+1]: %f\n", bpcov[sno][end+1]);


	if(neutral) { // if neutral (i.e. stranded) gets added to bpcov[1] here too => bpcov[1]=all coverage
		bpcov[1][start+1]+=v;
		bpcov[1][end+1]-=v;
	}
}


void add_read_to_cov_APPLY_UNISPG(GList<CReadAln>& rd,int n,GVec<float> *bpcov,int refstart, int refend) {

	fprintf(stderr, ">>> In add_read_to_cov, (int)rd[n]->strand: %d \n", (int)rd[n]->strand);

	int sno=(int)rd[n]->strand+1; // 0(-),1(.),2(+)

	fprintf(stderr, ">>> In add_read_to_cov, sno: %d \n", sno);
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
			fprintf(stderr, ">>> nsegs: %d \n", nsegs);

			int counter = 0;
			while(nsegs) {
				fprintf(stderr, ">>> loop %d ~~~\n", counter);
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
				fprintf(stderr, "\t>>> In add_read_to_cov, sno: %d \n", sno);
				fprintf(stderr, "\t>>> In add_read_to_cov, snop: %d \n", snop);
				if(sno!=snop) {
					if(sno==1) strand=snop;
					else if(snop!=1) strand=1;
				}
				fprintf(stderr, "\t>> Iterating pair_idx!!!\n");
				fprintf(stderr, "\t>>>>> start: %d;  refstart-refend: %d-%d;  start-refstart: %d\n", start, refstart, refend, start-refstart);
				fprintf(stderr, "\t>>>>> end: %d;  refstart-refend: %d-%d;  end-refstart: %d\n", end, refstart, refend, end-refstart);
				cov_edge_add_APPLY_UNISPG(bpcov,strand,start-refstart,end-refstart+1,rd[n]->pair_count[i]);
			}
		}
	}
	if(single_count>epsilon) {
		fprintf(stderr, ">> single_count>epsilon\n");
		for(int i=0;i<rd[n]->segs.Count();i++) {

			if (int(rd[n]->segs[i].end) < refend) continue;				

			int start = int(rd[n]->segs[i].start);
			int end = int(rd[n]->segs[i].end);
			if (start-refstart <= 0) {
				fprintf(stderr, "!!!! Warning!!!\n");
				start = refstart;
			}
			if (end >= refend) {
				end = refend;
			}
			fprintf(stderr, ">>>>> start: %u;  refstart: %d;  start-refstart: %d\n", rd[n]->segs[i].start, refstart, rd[n]->segs[i].start-refstart);
			fprintf(stderr, ">>>>> end: %u;  refstart: %d;  end-refstart: %d\n", rd[n]->segs[i].end, refstart, rd[n]->segs[i].end-refstart);
			cov_edge_add_APPLY_UNISPG(bpcov, sno, start-refstart, end-refstart+1, single_count);
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

	if((longreads||mixedMode) && bdata->keepguides.Count()) { // there are guides to consider with longreads -> I might need to adjust the splice sites
		GPVec<GffObj>& guides = bdata->keepguides;
		GList<CJunction> gjunc(true, true, true);
		for(int g=0;g<guides.Count();g++) {
			char s=0; // unknown strand
			if(guides[g]->strand=='+') s=1; // guide on positive strand
			else if(guides[g]->strand=='-') s=-1; // guide on negative strand
			for(int i=1;i<guides[g]->exons.Count();i++) { //add_junction((int)guides[g]->exons[i-1]->end,(int)guides[g]->exons[i]->start,gjunc,s);

				int gidx=gjunc.AddedIfNew(new CJunction(guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s));
				if(gidx>-1) { // item is new
					gjunc[gidx]->guide_match=1;
					if(mixedMode) { // this should favor junctions already covered by short reads
						CJunction jn(guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s);
						int oidx=-1;
						if (!junction.Found(&jn, oidx) || junction[oidx]->nm>=junction[oidx]->nreads) {
							gjunc[gidx]->guide_match=0;
							break;
						}
					}
				}
				fprintf(stderr,"Add guide junction=%d-%d:%d from reference %s\n",guides[g]->exons[i-1]->end,guides[g]->exons[i]->start,s,guides[g]->getID());
			}
		}

		if(gjunc.Count()) {

			//for(int i=0;i<gjunc.Count();i++) fprintf(stderr,"Guide junction=%d-%d:%d\n",gjunc[i]->start,gjunc[i]->end,gjunc[i]->strand);

			//int ssdist=1/ERROR_PERC;
			//int ssdist=longintronanchor;

			GVec<int> smodjunc; // keeps a list of all modified start junctions' indices (might keep the same junction twice)
			GVec<int> emodjunc; // keeps a list of all modified end junctions' indices (might keep the same junction twice)
			GList<CJunction> ejunction(junction);
			ejunction.setFreeItem(false);
			if(ejunction.Count()) ejunction.setSorted(juncCmpEnd);
			GList<CJunction> egjunc(gjunc);
			egjunc.setFreeItem(false);
			egjunc.setSorted(juncCmpEnd);
			int s=0;
			int e=0;
			for(int i=1;i<junction.Count();i++) {

				if(junction[i]->nm>=junction[i]->nreads){ // for all junctions -> try to see if I can correct them

					//fprintf(stderr,"check junction[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,junction[i]->start,junction[i]->end,junction[i]->strand,junction[i]->rightsupport,junction[i]->nm,junction[i]->nreads);

					// check start junction
					while(s<gjunc.Count() && gjunc[s]->start+sserror<junction[i]->start) s++;
					int k=s;
					int c=-1;
					int dist=1+sserror;
					while(k<gjunc.Count() && gjunc[k]->start<=junction[i]->start+sserror) {
						if(!junction[i]->strand || gjunc[k]->strand==junction[i]->strand) {
							if(gjunc[k]->start==junction[i]->start && gjunc[k]->guide_match) { // perfect match --> no need to change anything
								c=-1;
								break;
							}
							int d=dist;
							if(c<0 || gjunc[c]->guide_match==gjunc[k]->guide_match) d=abs((int)gjunc[k]->start-(int)junction[i]->start);
							if(d<dist || (!gjunc[c]->guide_match && gjunc[k]->guide_match)) {
								dist=d;
								c=k;
								smodjunc.Add(i);
							}
						}
						k++;
					}
					if(c>=0) {
						//fprintf(stderr,"...correct start of junction[%d] to %d\n",i,gjunc[c]->start);
						junction[i]->start=gjunc[c]->start;
						junction[i]->strand=gjunc[c]->strand;
						while(c<gjunc.Count() && gjunc[c]->start==junction[i]->start) {
							if(junction[i]->end==gjunc[c]->end && junction[i]->strand==gjunc[c]->strand) {
								junction[i]->guide_match=true;
								break;
							}
							c++;
						}
					}
				}
				if(ejunction[i]->nm>=ejunction[i]->nreads){ // for all junctions -> try to see if I can correct them

					//fprintf(stderr,"check ejunction[%d]:%d-%d:%d rightsupport=%f nm=%f nreads=%f\n",i,ejunction[i]->start,ejunction[i]->end,ejunction[i]->strand,ejunction[i]->rightsupport,ejunction[i]->nm,ejunction[i]->nreads);

					// check end junction
					while(e<egjunc.Count() && egjunc[e]->end+sserror<ejunction[i]->end) e++;
					int k=e;
					int c=-1;
					int dist=1+sserror;
					while(k<egjunc.Count() && egjunc[k]->end<=ejunction[i]->end+sserror) {
						if(!ejunction[i]->strand || egjunc[k]->strand==ejunction[i]->strand) {
							if(egjunc[k]->end==ejunction[i]->end  && egjunc[k]->guide_match) { // perfect match --> no need to change anything
								c=-1;
								break;
							}

							int d=dist;
							if(c<0 || egjunc[c]->guide_match==egjunc[k]->guide_match) d=abs((int)egjunc[k]->end-(int)ejunction[i]->end);
							if(d<dist || (!egjunc[c]->guide_match && egjunc[k]->guide_match)) {
								dist=d;
								c=k;
								emodjunc.Add(i);
							}
						}
						k++;
					}
					if(c>=0) {
						//fprintf(stderr,"...correct end of ejunction[%d] to %d\n",i,egjunc[c]->end);
						ejunction[i]->end=egjunc[c]->end;
						ejunction[i]->strand=egjunc[c]->strand;
						while(c<egjunc.Count() && egjunc[c]->end==ejunction[i]->end) {
							if(ejunction[i]->start==egjunc[c]->start && ejunction[i]->strand==egjunc[c]->strand) {
								ejunction[i]->guide_match=true; break;
							}
							c++;
						}
					}
					if(s==gjunc.Count() && e==egjunc.Count()) break;
				}
			}

			if(smodjunc.Count()||emodjunc.Count()) {
				modified=true;
				for(int i=0;i<smodjunc.Count();i++) {
					int j=smodjunc[i]; // junction that I modified --> try to see if there are other equal junctions
					//sprintf(sbuf, "%p", junction[j]);
					const CJunction* jp=jhash[junction[j]];
					if(!jp) { // did not process junction before
						//fprintf(stderr,"smodified junction[%d]:%d-%d:%d %p nreads=%f\n",j,junction[j]->start,junction[j]->end,junction[j]->strand,junction[j],junction[j]->nreads);
						s=j-1;
						GVec<int> equal;
						while(s>=0 && junction[s]->start==junction[j]->start) {
							if(junction[s]->end==junction[j]->end && junction[s]->strand==junction[j]->strand) equal.Add(s);
							s--;
						}
						s=j+1;
						while(s<junction.Count() && junction[s]->start==junction[j]->start) {
							if(junction[s]->end==junction[j]->end && junction[s]->strand==junction[j]->strand) equal.Add(s);
							s++;
						}
						if(equal.Count()) { // junction j is equal to other junctions
							for(s=0;s<equal.Count();s++) {
								//sprintf(sbuf, "%p", junction[equal[s]]);
								//jhash.Add(sbuf,junction[j]);
								CJunction* jct=junction[equal[s]];
								//fprintf(stderr,"...equal to junction[%d]:%d-%d:%d %p nreads=%f\n",equal[s],jct->start,jct->end,jct->strand,jct,jct->nreads);
								jhash.Add(jct, junction[j]);
								if(mixedMode) { // I can trust the reads coming from mixed data but not so much from long reads
									junction[j]->nreads+=jct->nreads;
									junction[j]->nm+=jct->nm;
									jct->nreads=0;
								}
								jct->strand=0;
								jct->guide_match=false;
							}
						}
					}
				}
				for(int i=0;i<emodjunc.Count();i++) {
					int j=emodjunc[i];
					//sprintf(sbuf, "%p", ejunction[j]);
					CJunction* jp=jhash[ejunction[j]];
					if(!jp) { // did not process junction before
						//fprintf(stderr,"emodified ejunction[%d]:%d-%d:%d %p nreads=%f\n",j,ejunction[j]->start,ejunction[j]->end,ejunction[j]->strand,ejunction[j],ejunction[j]->nreads);
						s=j-1;
						GVec<int> equal;
						while(s>=0 && ejunction[s]->end==ejunction[j]->end) {
							if(ejunction[s]->start==ejunction[j]->start && ejunction[s]->strand==ejunction[j]->strand) equal.Add(s);
							s--;
						}
						s=j+1;
						while(s<ejunction.Count() && ejunction[s]->end==ejunction[j]->end) {
							if(ejunction[s]->start==ejunction[j]->start && ejunction[s]->strand==ejunction[j]->strand) equal.Add(s);
							s++;
						}
						if(equal.Count()) { // junction j is equal to other junctions
							for(s=0;s<equal.Count();s++) {
								//sprintf(sbuf, "%p", ejunction[equal[s]]);
								//jhash.Add(sbuf,ejunction[j]);
								//ejunction[equal[s]]->strand=0;
								//ejunction[equal[s]]->guide_match=false;
								CJunction* ej=ejunction[equal[s]];
								jhash.Add(ej, ejunction[j]);
								//fprintf(stderr,"...equal to ejunction[%d]:%d-%d:%d %p nreads=%f\n",equal[s],ej->start,ej->end,ej->strand,ej,ej->nreads);
								if(mixedMode) { // I can trust the reads coming from mixed data but not so much from long reads
									ejunction[j]->nreads+=ej->nreads;
									ejunction[j]->nm+=ej->nm;
									ej->nreads=0;
								}
								ej->strand=0;
								ej->guide_match=false;
							}
						}
					}
				}
			}
			junction.Sort();
		}
		gjunc.Clear();
	}

	for(int s=0;s<3;s++) bpcov[s].Resize(refend-refstart+3);

	GVec<int> unstranded; // remembers unstranded reads

	/*****************************
	 * Iterating all reads
	 *****************************/
	for(int n=0;n<readlist.Count();n++) {
		CReadAln & rd=*(readlist[n]);

		fprintf(stderr, ">> rd[n]->pair_idx.Count(): %d\n", rd.pair_idx.Count());


		float rdcount=rd.read_count;

		int nex=rd.segs.Count();

		if(!rd.unitig) add_read_to_cov_APPLY_UNISPG(readlist,n,bpcov,refstart, refend);
		else if(rdcount>1) rdcount=1;

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
				/* DEL AWARE*/
				else if(!rd.juncs[i-1]->strand && (longreads || mixedMode) && (rd.segs[i-1].end!=rd.juncs[i-1]->start || rd.segs[i].start!=rd.juncs[i-1]->end)){ // see if I need to adjust read start/ends due to junction having deletions around it

					CJunction *nj=NULL;
					CJunction jn(rd.juncs[i-1]->start, rd.juncs[i-1]->end, rd.strand); // if the deletion makes sense then keep it
					int oidx=-1;
					if (junction.Found(&jn, oidx)) {
						nj=junction.Get(oidx);
						nj->nreads+=rd.juncs[i-1]->nreads;
						nj->mm+=rd.juncs[i-1]->mm;
						nj->nm+=rd.juncs[i-1]->nm;
						rd.juncs[i-1]=nj;
						// adjust start/end of read
						rd.segs[i-1].end=rd.juncs[i-1]->start;
						rd.segs[i].start=rd.juncs[i-1]->end;
					}
					else { // restore the correct junction
						CJunction jn(rd.segs[i-1].end, rd.segs[i].start, rd.strand);
						if (junction.Found(&jn, oidx)) {
							nj=junction.Get(oidx);
							nj->nreads+=rd.juncs[i-1]->nreads;
							nj->mm+=rd.juncs[i-1]->mm;
							nj->nm+=rd.juncs[i-1]->nm;
							rd.juncs[i-1]=nj;
						}
						else { // new junction
							rd.juncs[i-1]->strand=rd.strand;
							rd.juncs[i-1]->start=rd.segs[i-1].end;
							rd.juncs[i-1]->end=rd.segs[i].start;
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
			if(rd.juncs[i-1]->len()>longintron || (!longreads && rd.juncs[i-1]->nreads-rd.juncs[i-1]->nm<junctionthr && !rd.juncs[i-1]->mm)) {
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

		if((mixedMode || longreads) && !rd.strand && nex>1) unstranded.Add(n);

		//if(nex>1) fprintf(stderr,"\n");
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

	if(longreads || mixedMode) {
		for(int n=0;n<unstranded.Count();n++) {
			CReadAln & rd=*(readlist[unstranded[n]]);
			if(rd.start>=(uint)refstart && (int)(rd.end-(uint)refstart)<bpcov[1].Count()) {
				int nex=rd.segs.Count();
				float pos=0;
				float neg=0;
				if(longreads) {
					pos=get_cov(2,rd.segs[0].start-refstart,rd.segs[0].end-refstart,bpcov);
					neg=get_cov(0,rd.segs[0].start-refstart,rd.segs[0].end-refstart,bpcov);
				}
				bool posjunc=false;
				bool negjunc=false;
				for(int i=1;i<nex;i++) { // some junctions don't get an assigned strand by buggy aligners
					if(longreads) {
						pos+=get_cov(2,rd.segs[i].start-refstart,rd.segs[i].end-refstart,bpcov);
						neg+=get_cov(0,rd.segs[i].start-refstart,rd.segs[i].end-refstart,bpcov);
					}
					if(rd.juncs[i-1]->strand>0) posjunc=true;
					else if(rd.juncs[i-1]->strand<0) negjunc=true;
					else {
						int oidx=-1;
						CJunction jn(rd.juncs[i-1]->start, rd.juncs[i-1]->end, 1);
						if (junction.Found(&jn, oidx)) {
							posjunc=true;
						}
						jn.strand=-1;
						if (junction.Found(&jn, oidx)) {
							negjunc=true;
						}
					}
				}
				if(posjunc && !negjunc) {
					rd.strand=1;
				}
				else if(negjunc && !posjunc) {
					rd.strand=-1;
				}
				else {
					if(neg<1) neg=0;
					if(pos<1) pos=0;
					if(neg>pos) {rd.strand=-1;}
					else if(pos>neg) { rd.strand=1;}
				}
				//fprintf(stderr,"read strand is:%d for pos=%f neg=%f\n",rd.strand,pos,neg);
				if(rd.strand) {
					for(int i=1;i<nex;i++)
						if(!rd.juncs[i-1]->strand) {
							// first check to see if the junction already exists
							int oidx=-1;
							CJunction *nj=NULL;
							CJunction jn(rd.juncs[i-1]->start, rd.juncs[i-1]->end, rd.strand);
							if (junction.Found(&jn, oidx)) {
								nj=junction.Get(oidx);
								//fprintf(stderr,"found strand junction at %p",nj);
								nj->nreads+=rd.juncs[i-1]->nreads;
								nj->nreads_good+=rd.juncs[i-1]->nreads_good;
								nj->nm+=rd.juncs[i-1]->nm;
								nj->mm+=rd.juncs[i-1]->mm;
								nj->leftsupport+=rd.juncs[i-1]->leftsupport;
								nj->rightsupport+=rd.juncs[i-1]->rightsupport;
								rd.juncs[i-1]=nj;
							}
							else rd.juncs[i-1]->strand = rd.strand;
					}
				}
			}
		}
	}
}


CTransfrag *findtrf_in_treepat_APPLY_UNISPG(int gno,GIntHash<int>& gpos,GVec<int>& node,GBitVec& pattern,CTreePat *tr2no) { // doesn't work for patterns including source node
	fprintf(stderr, ">> Inside 'findtrf_in_treepat_APPLY_UNISPG'\n");
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