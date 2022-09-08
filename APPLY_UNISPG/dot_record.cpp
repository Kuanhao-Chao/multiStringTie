#include "dot_record.h"
#include "GBitVec.h"
#include "unispg.h"
#include <float.h>

#define GMEMTRACE 1  //debugging memory allocation
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

bool DOTInputFile::start(const char* fn) {

	fprintf(stderr, "Inside start %s\n", fn);
    this -> file = fn;
    this -> tmpfile = fn;
	//stringtie multi-BAM input
	// for (int i=0;i<bamfiles.Count();++i) {
	this -> reader = new DOTReader(this -> file);
    bool is_dot_open;

	is_dot_open = this -> reader -> dotopen(this -> file);
	fprintf(stderr, "is_dot_open: %d\n", is_dot_open);
	// if (is_dot_open) {
	// 	this -> rec = this -> reader->next();
	// }
	return is_dot_open;
    // UnispgGp_APPLY* brec=this -> reader->next();
}


UnispgGp_APPLY* DOTInputFile::next() {
	//must free old current record first
	// delete rec;
	rec=NULL;
	rec = reader->next();
	// this -> updateUnispgGp_APPLYGp();
	// uni_splice_graphs
	return rec;
    // if (recs.Count()>0) {
    // 	crec=recs.Pop();//lowest coordinate
    	// UnispgGp_APPLY* rnext = this->reader->next();
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
    // unlink(tmpfile.chars());
//  }
}

