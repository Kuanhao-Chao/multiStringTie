#include "tree_pattern_A.h"

CTreePat *construct_treepat_APPLY_UNISPG(int node_num, GIntHash<int>& gpos, GPVec<CTransfrag>& transfrag) {
	// fprintf(stderr, ">> Inside 'construct_treepat_APPLY_UNISPG' (%d) \n", node_num);
	// create root CTreePat first
	CTreePat *root=new CTreePat(0,node_num-1); // if links from source to nodes are desired source==1 and all nodes are treated as +1

	// now construct all child CTreePat's
	// fprintf(stderr,"There are %d transfrags\n",transfrag.Count());
	for(int t=0;t<transfrag.Count();t++)
		if(transfrag[t]->nodes[0]){ // don't include transfrags from source -> not needed
			CTreePat *tree=root;
			int m=0; // previous node in pattern that was set in pattern
			for(int n=1;n<node_num;n++){
				if(transfrag[t]->pattern[n]) {
					CTreePat *child;
					if(m) { // there is a node m that was seen before
						int *pos=gpos[edge(m,n,node_num)];
						if(pos && transfrag[t]->pattern[*pos]) // there is an edge between m and n
							child=tree->settree(node_num-1-m+n-m-1,n,2*(node_num-n-1));
						else child=tree->settree(n-m-1,n,2*(node_num-n-1));
					}
					else // this is the root tree
						child=tree->settree(n-1,n,2*(node_num-n-1));
					tree=child;
					m=n;
				}
			}
			tree->tr=transfrag[t];
		}

	return(root);
}
