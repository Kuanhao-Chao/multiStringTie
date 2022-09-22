#ifndef __STATIC_FUNC_H__
#define __STATIC_FUNC_H__

static inline int edge(int min, int max, int gno) {
	//return((gno-1)*min-min*(min-1)/2+max-min); // this should be changed if source to node edges are also stored
	return((gno-1)*(min+1)-min*(min-1)/2+max-min); // this includes source to node edges
}

#endif