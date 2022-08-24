#include "parse_reads.h"

void get_fragment_pattern(GList<CReadAln>& readlist, int n, int np, GPVec<UnispgGp>** graphs_vec) {
    fprintf(stderr, "get_fragment_pattern is called. \n");

	/*****************************
	 * Step 1: Generate the read pattern !!!
	 *  both reads are unstranded
	 *****************************/
	for(int s=0;s<2;s++) {

        /*****************************
         * Step 1-1: After the `get_read_pattern`, 
         *   (1) confirm each read belongs to which nodes
         *   (2) the node coverage gonna be updated
         *****************************/
    	// rgno: It stores the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> rgno;
	    // rnode: It stores the node indices in the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> *rnode=new GVec<int>[readlist[n]->segs.Count()];

		if(readlist[n]->nh) {
			get_read_pattern(s, readlist, n, rgno, rnode, graphs_vec);
        }



    	// pgno: It stores the paired bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> pgno;
	    // pnode: It stores the paired node indices in the bundle<->graph indices in a BundleData that a read belongs to.
		GVec<int> *pnode=NULL;

        if(np>-1 && readlist[np]->nh) {
			pnode=new GVec<int>[readlist[np]->segs.Count()];
			get_read_pattern(s, readlist, n, rgno, rnode, graphs_vec);
		}


        /*****************************
         * Step 1-2: actually create a pattern for a read
         *****************************/
        // Read pattern
        GBitVec rpat;
        // Paired pattern
        GBitVec ppat;
        int usedp=0;

        // set read pattern
        rpat.clear();
        // rpat.resize(graphno[s][rgno[r]]+edgeno[s][rgno[r]]);
        int i=0;
        int clear_rnode=0;

    }
}


/*****************************
 * (1) confirm each read belongs to which nodes
 * (2) the node coverage gonna be updated
 *****************************/
void get_read_pattern(int s, GList<CReadAln>& readlist, int n,GVec<int> rgno, GVec<int> *rnode, GPVec<UnispgGp>** graphs_vec) {

	// rgno: It stores the bundle<->graph indices in a BundleData that a read belongs to.
	// rnode: It stores the node indices in the bundle<->graph indices in a BundleData that a read belongs to.

	int lastgnode=-1;
	int lastngraph=-1;
	int ncoord=readlist[n]->segs.Count();
	int k=0; // need to keep track of coordinates already added to coverages of graphnodes
	int kmer=KMER-1; //f1

	GIntHash<bool> hashnode;

    // Here, I need to link a read to the node/nodes that it belongs to

}