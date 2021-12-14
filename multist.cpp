#include "multist.h"
#include "GBitVec.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

bool DOTInputFile::start(const char* fn) {
    this -> file = fn;
    this -> tmpfile = fn;
	//stringtie multi-BAM input
	// for (int i=0;i<bamfiles.Count();++i) {
	this -> reader = new DOTReader(this -> file);
    bool is_dot_open;
	is_dot_open = this -> reader -> dotopen(this -> file);
	// if (is_dot_open) {
	// 	this -> rec = this -> reader->next();
	// }
	return is_dot_open;

    // DOTRecord* brec=this -> reader->next();
}


DOTRecord* DOTInputFile::next() {
	//must free old current record first
	delete this -> rec;
	this -> rec=NULL;
	this -> rec = this->reader->next();
	return this -> rec;
    // if (recs.Count()>0) {
    // 	crec=recs.Pop();//lowest coordinate
    	// DOTRecord* rnext = this->reader->next();
    	// if (rnext)
    	// 	recs.Add(new DOTInputRecord(rnext, crec->fidx));
    	// return crec->brec;
		// return rnext;
    // }
    // else return NULL;
}


void DOTInputFile::stop() {
	delete rec;
	rec=NULL;
 	this->reader -> dotclose();
//  if (!keepTempFiles) {
    unlink(tmpfile.chars());
//  }
}

