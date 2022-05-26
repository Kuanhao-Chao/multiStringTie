#include "rlink.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif


#define BSIZE 10000 // bundle size

static GStr _id("", 256); //to prevent repeated reallocation for each parsed read

bool exonmatch(GVec<GSeg> &prevexons, GVec<GSeg> &exons) {
	if(prevexons.Count() != exons.Count()) return false;
	for(int i=0;i<exons.Count();i++) {
		if(prevexons[i].end!=exons[i].end || prevexons[i].start!=exons[i].start) return false;
	}
	return true;
}


bool mismatch_anchor(CReadAln *rd,char *mdstr,int refstart, bam1_t *b) {
	if(mdstr==NULL) return false;
	//--make a copy of the string, in case the original is a const string
	// (because the parseUInt() function modifies the string temporarily
	char* mdstring=Gstrdup(mdstr);
	char *p=mdstring;

	int i=0;
	int parsedlen=refstart;
	int rdlen=0;

	while (*p!='\0') {
		unsigned int num_matches=0;
		if (*p>='0' && *p<='9') {
			parseUInt(p,num_matches);
			//if (num_matches>0) GMessage("%d matching bases\n", num_matches);
			parsedlen+=num_matches;
			continue;
		}
		if (*p=='^') { //deletion --> found a problem with deletion
			//GDynArray<char> deletion; //deletion string accumulates here (if needed)
			int del_length=0;//tracking deletion length
			char delbase=*(++p);
			while (delbase>='A' && delbase<='Z') {
				//deletion.Add(delbase);
				del_length++;
				delbase=*(++p);
			}
			while(i<rd->segs.Count() && rdlen+(int)rd->segs[i].len()<parsedlen) {
				rdlen+=rd->segs[i].len();
				i++;
			}
			if(i==rd->segs.Count()) break;

			if((i && parsedlen-rd->segs[i].start<junctionsupport) || (i<rd->segs.Count()-1 && rd->segs[i].end+1-parsedlen-del_length<junctionsupport)) {
				GFREE(mdstring);
				return true;
			}
			parsedlen+=del_length;

			/*GMessage("%d base(s) deletion [", del_length);
			for (uint i=0;i<deletion.Count();++i) GMessage("%c",deletion[i]);
			GMessage("]\n");*/
			continue;
		}
		if (*p>='A' && *p<='Z') {
			//GMessage("base mismatch [%c]\n",*p);
			while(i<rd->segs.Count() && rdlen+(int)rd->segs[i].len()<parsedlen) {
				rdlen+=rd->segs[i].len();
				i++;
			}
			if(i==rd->segs.Count()) break;
			if((i && parsedlen-rd->segs[i].start<junctionsupport) || (i<rd->segs.Count()-1 && rd->segs[i].end-parsedlen<junctionsupport)) {
				GFREE(mdstring);
				return true;
			}
			parsedlen++;
		}
		p++;
	}

	GFREE(mdstring);

	uint32_t *cigar = bam_get_cigar(b);
	rdlen=0;
	parsedlen=0;
	i=0;

	for (uint j = 0; j < b->core.n_cigar; ++j) {
		int op = bam_cigar_op(cigar[j]);
		if (op == BAM_CMATCH || op==BAM_CEQUAL ||
				op == BAM_CDIFF || op == BAM_CDEL) {
			parsedlen += bam_cigar_oplen(cigar[j]);
		}
		else if(op == BAM_CINS) {
			while(i<rd->segs.Count() && rdlen+(int)rd->segs[i].len()<parsedlen) {
				rdlen+=rd->segs[i].len();
				i++;
			}
			if(i==rd->segs.Count()) break;
			if((i && parsedlen-rd->segs[i].start<junctionsupport) || (i<rd->segs.Count()-1 && rd->segs[i].end-parsedlen<junctionsupport)) {
				return true;
			}
		}
	}

	return false;
}

void processRead(int currentstart, int currentend, BundleData& bdata,
		 GHash<int>& hashread,  GReadAlnData& alndata) { // some false positives should be eliminated here in order to break the bundle

	// fprintf(stderr, ">> Inside processRead !! currentstart: %d, currentend: %d \n", currentstart, currentend);
	GSamRecord& brec=*(alndata.brec);			   // bam record
	GList<CReadAln>& readlist = bdata.readlist;    // list of reads gathered so far
	GList<CJunction>& junction = bdata.junction;   // junctions added so far
    char strand=alndata.strand;
    int nh=alndata.nh;
    int hi=alndata.hi;
	int readstart=brec.start;
	CReadAln* readaln=NULL;                        // readaln is initialized with NULL
	//bool covSaturated=false;                       // coverage is set to not saturated


	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Process read %s with strand=%d and exons:",brec.name(),strand);
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
		}
		fprintf(stderr,"\n");
	}
	*/

	double nm=(double)brec.tag_int("NM"); // read mismatch
	float unitig_cov=0;
	unitig_cov=brec.tag_float("YK");

	bool longr=false;

	bool match=false;  // true if current read matches a previous read
	int n=readlist.Count()-1;

	while(n>-1 && readlist[n]->start==brec.start) {
		if(strand==readlist[n]->strand && (readlist[n]->longread==longr) && (!isunitig || (unitig_cov>0) == readlist[n]->unitig)) {
			match=exonmatch(readlist[n]->segs,brec.exons);
		}
		//if(strand==readlist[n]->strand) match=exonmatch(readlist[n]->segs,brec.exons);
		if(match) break; // this way I make sure that I keep the n of the matching readlist
		n--;
	}

	if (bdata.end<currentend) {// I am not sure why this is done here?
		bdata.start=currentstart;
		bdata.end=currentend;
	}
	bdata.numreads++;                         // number of reads gets increased no matter what
	//bdata.wnumreads+=float(1)/nh;

	if (!match) { // if this is a new read I am seeing I need to set it up
		readaln=new CReadAln(strand, nh, brec.start, brec.end, alndata.tinfo);
		readaln->longread=longr;
		alndata.tinfo=NULL; //alndata.tinfo was passed to CReadAln
		for (int i=0;i<brec.exons.Count();i++) {
			readaln->len+=brec.exons[i].len();
			if(i) {
				if(!junction.Count()) { // always add null junction first
					CJunction *nullj=new CJunction(0, 0, 0);
					junction.Add(nullj);
				}
				int jstrand=strand;
				uint jstart=brec.exons[i-1].end;
				uint jend=brec.exons[i].start;

				//fprintf(stderr,"exon count=%d junctiondel count=%d exonend=%d exonstart=%d\n",brec.exons.Count(),brec.juncsdel.Count(),jstart,jend);

				CJunction* nj=junction.AddIfNew(new CJunction(jstart, jend, jstrand), true);
				if (alndata.juncs.Count())
					nj->guide_match=alndata.juncs[i-1]->guide_match;
				if (nj) {
					readaln->juncs.Add(nj);
				}
			}
			readaln->segs.Add(brec.exons[i]);
		}
		n=readlist.Add(readaln);  // reset n for the case there is no match
	}
	else { //redundant read alignment matching a previous alignment
		// keep shortest nh so that I can see for each particular read the multi-hit proportion
		if(nh<readlist[n]->nh) readlist[n]->nh=nh;
		/*
		//for mergeMode, we have to free the transcript info:
		if (alndata.tinfo!=NULL) {
			 delete alndata.tinfo;
			 alndata.tinfo=NULL;
		}
		*/
	}

	/*
	{ // DEBUG ONLY
		fprintf(stderr,"Add read %s with strand=%d and exons:",brec.name(),strand);
		for (int i=0;i<brec.exons.Count();i++) {
			fprintf(stderr," %d-%d", brec.exons[i].start, brec.exons[i].end);
		}
		fprintf(stderr,"\n");
		//fprintf(stderr,"Read %s is at n=%d with unitig_cov=%f and strand=%d\n",brec.name(),n,unitig_cov,strand);
	}
	*/

	if((int)brec.end>currentend) {
		currentend=brec.end;
	  	bdata.end=currentend;
	}

	float rdcount=(float)brec.tag_int("YC"); // alignment count
	if(!rdcount) rdcount=1;
	if(unitig_cov) {
		rdcount=unitig_cov;
		if(isunitig) readlist[n]->unitig=true; // treat unitig reads differently when -U is set
	}

	if(!nomulti) rdcount/=nh;
	readlist[n]->read_count+=rdcount; // increase single count just in case I don't see the pair

	// store the mismatch count per junction so that I can eliminate it later
	if(!nm) {
		nm=(double)brec.tag_int("nM"); // paired mismatch : big problem with STAR alignments
		if(brec.isPaired()) nm/=2;
	}
	if(brec.clipL) nm++;
	if(brec.clipR) nm++;

	//nm+=brec.clipL; Note: this clippings were to aggressive
	//nm+=brec.clipR;

	if(readlist[n]->juncs.Count()) {
		bool mismatch=false;
		if(readlist[n]->longread) mismatch=true;
		else if(nm/readlist[n]->len>mismatchfrac) mismatch=true;
		else if(nm && readlist[n]->juncs.Count()) {
			if(brec.clipL && readlist[n]->segs[0].len()<junctionsupport+brec.clipL) mismatch=true; // penalize mismatch that's too close to ss
			else if(brec.clipR && readlist[n]->segs.Last().len()<junctionsupport+brec.clipR) mismatch=true;
			else if(mismatch_anchor(readlist[n],brec.tag_str("MD"),currentstart,brec.get_b())) mismatch=true; // this line was not initially present in vs1 or vs3 but I noticed it doesn't do any difference in real data, so far it only helped with the SR in simulation -> I might want to take it out
		}

		for(int i=0;i<readlist[n]->juncs.Count();i++) { // if read is PacBio I might want to increase the mismatch fraction, although the nm only gets used for longintrons
			if(mismatch || nh>2) readlist[n]->juncs[i]->nm+=rdcount;
			if(readlist[n]->segs[i].len()>longintronanchor && readlist[n]->segs[i+1].len()>longintronanchor)
				readlist[n]->juncs[i]->mm+=rdcount;
			//if(nh>2) readlist[n]->juncs[i]->mm+=rdcount; // vs3; vs2 only has nh>1
			readlist[n]->juncs[i]->nreads+=rdcount;
		}
	}


	// now set up the pairing
	if (brec.refId()==brec.mate_refId()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
	//if (brec.refId()==brec.mate_refId() && brec.isProperlyPaired()) {  //only consider mate pairing data if mates are on the same chromosome/contig and are properly paired
	//if (brec.isProperlyPaired()) {  //only consider mate pairing data if mates  are properly paired
		int pairstart=brec.mate_start();
		if (currentstart<=pairstart) { // if pairstart is in a previous bundle I don't care about it
			//GStr readname();
			//GStr id(brec.name(), 16); // init id with readname
			_id.assign(brec.name()); //assign can be forced to prevent shrinking of the string
			if(pairstart<=readstart) { // if I've seen the pair already <- I might not have seen it yet because the pair starts at the same place
				_id+='-';_id+=pairstart;
				_id+=".=";_id+=hi; // (!) this suffix actually speeds up the hash by improving distribution!
				const int* np=hashread[_id.chars()];
				if(np) { // the pair was stored --> why wouldn't it be? : only in the case that the pair starts at the same position
					if(readlist[*np]->nh>nh && !nomulti) rdcount=float(1)/readlist[*np]->nh;
					bool notfound=true;
					for(int i=0;i<readlist[*np]->pair_idx.Count();i++)
						if(readlist[*np]->pair_idx[i]==n) {
							readlist[*np]->pair_count[i]+=rdcount;
							notfound=false;
							break;
						}
					if(notfound) { // I didn't see the pairing before
						readlist[*np]->pair_idx.Add(n);
						readlist[*np]->pair_count.Add(rdcount);
					}

					notfound=true;
					for(int i=0;i<readlist[n]->pair_idx.Count();i++)
						if(readlist[n]->pair_idx[i]==*np) {
							readlist[n]->pair_count[i]+=rdcount;
							notfound=false;
							break;
						}
					if(notfound) { // I didn't see the pairing before
						int i=*np;
						readlist[n]->pair_idx.Add(i);
						readlist[n]->pair_count.Add(rdcount);
					}
					hashread.Remove(_id.chars());
				}
			}
			else { // I might still see the pair in the future
				_id+='-';_id+=readstart; // this is the correct way
				_id+=".=";_id+=hi;
				hashread.Add(_id.chars(), n);
			}
		}
	} //<-- if mate is mapped on the same chromosome
	// fprintf(stderr, "Process read %d - %d \n", alndata.brec->start,alndata.brec->end);
	// fprintf(stderr, "Current start & end %d - %d \n", currentstart, currentend);
}


// get_fragment_pattern(bundle->readlist,n,np,bundle->readlist[n]->pair_count[j],readgroup,merge,group2bundle,bundle2graph,graphno,edgeno,gpos,no2gnode,transfrag,tr2no,group);

void get_fragment_pattern(GList<CReadAln>& readlist,int n, int np,float readcov,GVec<int> *readgroup,GVec<int>& merge, GVec<int> *group2bundle,
		GVec<CGraphinfo> **bundle2graph,GVec<int> *graphno,GVec<int> *edgeno, GIntHash<int> **gpos,GPVec<CGraphnode> **no2gnode,
		GPVec<CTransfrag> **transfrag,CTreePat ***tr2no,GPVec<CGroup> &group) {
	uint rstart=readlist[n]->start; // this only works for unpaired long reads -> I will have to take into account the pair if I want to do it for all reads
	uint rend=readlist[n]->end;
	if(np>-1 && readlist[np]->end>rend) rend=readlist[np]->end;

	float rprop[2]={0,0}; // by default read does not belong to any strand
	// compute proportions of unstranded read associated to strands

	if(readlist[n]->nh && !readlist[n]->strand && np>-1 && readlist[np]->nh && !readlist[np]->strand) { // both reads are unstranded
		
	} else {

	}
}