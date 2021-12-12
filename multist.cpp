#include "multist.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

void DOTInputFile::start() {
	GStr dotfile;
    dotfile = this -> file;

	//stringtie multi-BAM input
	// for (int i=0;i<bamfiles.Count();++i) {
    DOTReader* dotreader=new DOTReader(dotfile);
    DOTRecord* brec=dotreader->next();
}


DOTRecord* DOTInputFile::next() {
	//must free old current record first
	delete crec;
	crec=NULL;
    // if (recs.Count()>0) {
    // 	crec=recs.Pop();//lowest coordinate
    	DOTRecord* rnext = this->reader->next();
    	// if (rnext)
    	// 	recs.Add(new DOTInputRecord(rnext, crec->fidx));
    	// return crec->brec;
		return rnext;
    // }
    // else return NULL;
}


void DOTInputFile::stop() {
 	this->reader -> dotclose();

//  if (!keepTempFiles) {
    unlink(tmpfile.chars());
//  }
}

