#ifndef __T_RECORD_H__
#define __T_RECORD_H__

#pragma once

#include "global_params.h"

// #include "GStr.h"
// #include "GSam.h"
// #include "gff.h"
// #include "GList.hh"
// #include "rlink.h"
// extern GStr tmp_path;
// extern GStr cram_ref;
// extern bool keepTempFiles;
// extern bool debugMode; // "debug" or "D" tag.
// extern bool verbose; // "v" tag.

struct TInputRecord {
	GSamRecord* brec;
	int fidx; //index in files and readers
	bool operator<(TInputRecord& o) {
		 //decreasing location sort
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 //int refcmp=strcmp(r1.refName(),r2.refName());
		 int refcmp=r1.refId()-r2.refId();
		 if (refcmp==0) {
		 //higher coords first
			if (r1.start!=r2.start)
				 return (r1.start>r2.start);
			else {
				if (r1.end!=r2.end)
				   return (r1.end>r2.end);
				else if (fidx==o.fidx)
						return (strcmp(r1.name(), r2.name())>0);
					else return fidx>o.fidx;
			}
		 }
		 else { //use header order
			 return (refcmp>0);
		 }
	}
	bool operator==(TInputRecord& o) {
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 return ( strcmp(r1.refName(),r2.refName())==0 && r1.start==r2.start && r1.end==r2.end
				 && fidx==o.fidx && strcmp(r1.name(),r2.name())==0);
	}

	TInputRecord(GSamRecord* b=NULL, int i=0):brec(b),fidx(i) {}
	~TInputRecord() {
		delete brec;
	}
};

struct TInputFiles {
 protected:
	TInputRecord* crec;
	GStr convert2BAM(GStr& gtf, int idx);
 public:
	GPVec<GSamReader> readers;
	GVec<GStr> files; //same order
	GVec<GStr> tmpfiles; //all the temp files created by this
	
	// KH add
	GVec<GStr> bamfiles; //Processed file name

	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), readers(true), files(), tmpfiles(),
			recs(true, true, true) { 
				// fprintf(stderr, "TInputFiles initialization\n");
			}
	void Add(const char* fn);
	int count() { return files.Count(); }
	int start(); //open all files, load 1 record from each
	void start_fidx(int fidx);
	GSamRecord* next();
	void stop(); //
	void stop_fidx(int fidx);
};


#endif /* __T_RECORD_H__ */