#include "t_record.h"

void TInputFiles::Add(const char* fn) {
	   GStr sfn(fn);
		if (sfn!="-" && fileExists(fn)<2) {
			    GError("Error: input file %s cannot be found!\n",
			            fn);
		}

		files.Add(sfn);
	}


int gseqstat_cmpName(const pointer p1, const pointer p2) {
	return strcmp(((GSeqStat*)p1)->gseqname, ((GSeqStat*)p2)->gseqname);
}

GStr TInputFiles::convert2BAM(GStr& gtf, int idx) {
  GStr bamfname(tmp_path);
  bamfname.appendfmt("transcripts_s%04d",idx);
  GStr samhname(bamfname);
  bamfname+=".bam";
  samhname+=".sam";
  tmpfiles.Add(bamfname);
  tmpfiles.Add(samhname);
  FILE* samh=fopen(samhname.chars(), "w");
  if (samh==NULL) GError("Error creating file: %s\n",samhname.chars());
  fprintf(samh, "@HD\tVN:1.0\tSO:coordinate\n");
  //load GTF as sorted
  GffReader gfr(gtf.chars(), true, true); //transcript only, sorted by location
  gfr.setRefAlphaSorted(true); //make sure refseq IDs are sorted alphabetically
  gfr.showWarnings(debugMode || verbose);
  gfr.readAll(true, true, true); //keep attributes, merge close exons, no_exon_attributes
  if (gfr.gflst.Count()==0)
	  GError("Error: no transcripts were found in input file %s\n", gtf.chars());
  gfr.gseqStats.Sort(gseqstat_cmpName);
  for (int i=0;i<gfr.gseqStats.Count();++i) {
  	fprintf(samh, "@SQ\tSN:%s\tLN:%u\n", gfr.gseqStats[i]->gseqname,
  			gfr.gseqStats[i]->maxcoord+500);
  }
  fprintf(samh, "@CO\tfn:%s\n",gtf.chars());
  fclose(samh);
  GSamWriter bw(bamfname.chars(),samhname.chars());
  for (int i=0;i<gfr.gflst.Count();++i) {
	  GffObj& m = *gfr.gflst[i];
	  int t_id=bw.get_tid(m.getGSeqName());
	  if (t_id<0)
		   GError("Error getting header ID# for gseq %s (file: %s)\n",m.getGSeqName(),gtf.chars());
	  GDynArray<uint32_t> cigar;
	  for (int k=0;k<m.exons.Count();++k) {
		  if (k>0) {
			  cigar.Add((int(m.exons[k]->start-m.exons[k-1]->end-1) << BAM_CIGAR_SHIFT) | BAM_CREF_SKIP); // cigar+='N';
		  }
		  //cigar+=m.exons[k]->len();
		  //cigar+='M';
		  cigar.Add( (m.exons[k]->len() << BAM_CIGAR_SHIFT) | BAM_CMATCH);
	  }
	  GSamRecord brec(m.getID(), t_id, m.start, false, cigar);
	  if (m.strand=='-' || m.strand=='+') {
		   GStr tag("XS:A:");
		   tag+=m.strand;
	       brec.add_aux(tag.chars());
	  }
	  //GStr s("ZF:i:");
	  //s+=idx;
	  //brec.add_aux(s.chars());
	  char *av=m.getAttr("cov");
	  if (av!=NULL) {
		  GStr s("ZS:Z:",20);
		  s+=av;
		  s+='|';
		  av=m.getAttr("FPKM");
		  if (av) s+=av;
		  av=m.getAttr("TPM");
		  if (av) { s+='|';s+=av; }
		  brec.add_aux(s.chars());
	  }
	  bw.write(&brec);
  } //for each transcript
  return bamfname;
}


int TInputFiles::start() {
	// GVec<GStr> bamfiles;
	this->bamfiles=files;
	// //stringtie multi-BAM input
	// for (int i=0;i<bamfiles.Count();++i) {
	// 	GSamReader* bamreader=new GSamReader(bamfiles[i].chars(), cram_ref.is_empty() ? NULL : cram_ref.chars());
	// 	readers.Add(bamreader);
	// 	GSamRecord* brec=bamreader->next();
	// 	if (brec)
	// 	   {	
	// 		//    TInputRecord* test = new TInputRecord(brec, i);
	// 		   recs.Add(new TInputRecord(brec, i));
	// 		   	// fprintf(stderr, "test->fidx: %d\n", test->fidx);
	// 	   }
	// }
	return bamfiles.Count();
}


void TInputFiles::start_fidx(int fidx) {
	//stringtie multi-BAM input
	// for (int i=0;i<bamfiles.Count();++i) {
	GSamReader* bamreader=new GSamReader(this->bamfiles[fidx].chars(), cram_ref.is_empty() ? NULL : cram_ref.chars());
	readers.Add(bamreader);
	GSamRecord* brec=bamreader->next();
	if (brec) {	
		//    TInputRecord* test = new TInputRecord(brec, i);
		recs.Add(new TInputRecord(brec, fidx));
		// fprintf(stderr, "test->fidx: %d\n", test->fidx);
	}
	// }
	// return readers.Count();
}


GSamRecord* TInputFiles::next() {
	//must free old current record first
	delete crec;
	crec=NULL;
	// fprintf(stderr, "recs.Count(): %d\n", recs.Count());
    if (recs.Count()>0) {
    	crec=recs.Pop();//lowest coordinate
		// fprintf(stderr, "crec->fidx: %d\n", crec->fidx);
    	GSamRecord* rnext=readers[crec->fidx]->next();
		// fprintf(stderr, "rnext: %d - %d \n", rnext->start, rnext->end);
    	if (rnext)
    		recs.Add(new TInputRecord(rnext, crec->fidx));
    	crec->brec->uval=crec->fidx; //send file index
    	return crec->brec;
    }
    else return NULL;
}

// GSamRecord* TInputFiles::next_fidx(int fidx) {
// 	//must free old current record first
// 	delete crec;
// 	crec=NULL;
// 	// fprintf(stderr, "recs.Count(): %d\n", recs.Count());
//     if (recs.Count()>0) {
//     	crec=recs.Pop();//lowest coordinate
// 		// fprintf(stderr, "crec->fidx: %d\n", crec->fidx);
//     	GSamRecord* rnext=readers[crec->fidx]->next();
// 		// fprintf(stderr, "rnext: %d - %d \n", rnext->start, rnext->end);
//     	if (rnext)
//     		recs.Add(new TInputRecord(rnext, crec->fidx));
//     	crec->brec->uval=crec->fidx; //send file index
//     	return crec->brec;
//     }
//     else return NULL;
// }

void TInputFiles::stop() {
 for (int i=0;i<readers.Count();++i) {
	 readers[i]->bclose();
 }
 if (!keepTempFiles) {
	 for (int i=0;i<tmpfiles.Count();++i) {
		 unlink(tmpfiles[i].chars());
	 }
 }
}

void TInputFiles::stop_fidx(int fidx) {
 for (int i=0;i<readers.Count();++i) {
	 readers[i]->bclose();
 }
}