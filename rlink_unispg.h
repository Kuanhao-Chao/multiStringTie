#include "unispg.h"
/****************
 **  KH Adding 
****************/
extern FILE* uinigraph_out;
// extern bool universal_splice_graph;
extern GStr outfname_prefix;

extern FILE* node_cov_pos_bed;
extern FILE* edge_cov_pos_bed;
extern FILE* node_cov_neg_bed;
extern FILE* edge_cov_neg_bed;


// extern UnispgGp* unispg_gp;
extern GPVec<CGraphnode>* no2gnodeGp_unispg[2];
extern GVec<int> current_gidx;
extern int track_idx;

extern GVec<FILE*> node_cov_pos_bed_vec;
// GVec<GStr> nodecovposfname_vec; 
extern GVec<FILE*> edge_cov_pos_bed_vec;
// GVec<GStr> edgecovposfname_vec; 

extern GVec<FILE*> node_cov_neg_bed_vec;
// GVec<GStr> nodecovnegfname_vec; 
extern GVec<FILE*> edge_cov_neg_bed_vec;
// GVec<GStr> edgecovnegfname_vec; 

extern GVec<FILE*> node_cov_pos_bed_unispg_vec;
// GVec<GStr> nodecovposfname_unispg_vec; 
extern GVec<FILE*> edge_cov_pos_bed_unispg_vec;
// GVec<GStr> edgecovposfname_unispg_vec; 

extern GVec<FILE*> node_cov_neg_bed_unispg_vec;
// GVec<GStr> nodecovnegfname_unispg_vec; 
extern GVec<FILE*> edge_cov_neg_bed_unispg_vec;
// GVec<GStr> edgecovnegfname_unispg_vec; 
/****************
 **  END KH Adding 
****************/
// #include "visualization.h"

int infer_transcripts_unispg(BundleData* bundle, int fidx);