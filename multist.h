#ifndef __MULTIST_H__
#define __MULTIST_H__
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "time.h"
#include "tablemaker.h"
#include "GHashMap.hh"

// struct DOTInputRecord {
// 	GSamRecord* brec;
// 	int fidx; //index in files and readers
// 	bool operator<(TInputRecord& o) {
// 		 //decreasing location sort
// 		 GSamRecord& r1=*brec;
// 		 GSamRecord& r2=*(o.brec);
// 		 //int refcmp=strcmp(r1.refName(),r2.refName());
// 		 int refcmp=mergeMode ? strcmp(r1.refName(),r2.refName()) : r1.refId()-r2.refId();
// 		 if (refcmp==0) {
// 		 //higher coords first
// 			if (r1.start!=r2.start)
// 				 return (r1.start>r2.start);
// 			else {
// 				if (r1.end!=r2.end)
// 				   return (r1.end>r2.end);
// 				else if (fidx==o.fidx)
// 						return (strcmp(r1.name(), r2.name())>0);
// 					else return fidx>o.fidx;
// 			}
// 		 }
// 		 else { //use header order
// 			 return (refcmp>0);
// 		 }
// 	}
// 	bool operator==(TInputRecord& o) {
// 		 GSamRecord& r1=*brec;
// 		 GSamRecord& r2=*(o.brec);
// 		 return ( strcmp(r1.refName(),r2.refName())==0 && r1.start==r2.start && r1.end==r2.end
// 				 && fidx==o.fidx && strcmp(r1.name(),r2.name())==0);
// 	}

// 	TInputRecord(GSamRecord* b=NULL, int i=0):brec(b),fidx(i) {}
// 	~TInputRecord() {
// 		delete brec;
// 	}
// };

// struct DOTInputFile {
//  protected:
// 	TInputRecord* crec;
//  public:
// 	GPVec<GSamReader> readers;
// 	GVec<GStr> files; //same order
// 	GVec<GStr> tmpfiles; //all the temp files created by this
// 	GList<TInputRecord> recs; //next record for each
// 	TInputFiles():crec(NULL), readers(true), files(), tmpfiles(),
// 			recs(true, true, true) { }
// 	void Add(const char* fn);
// 	int count() { return files.Count(); }
// 	int start(); //open all files, load 1 record from each
// 	GSamRecord* next();
// 	void stop(); //
// };



#endif
