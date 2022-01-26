#ifndef __UNISPG_H__
#define __UNISPG_H__
#include "rlink.h"
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "time.h"
#include "tablemaker.h"
#include "GHashMap.hh"
#include <fstream>
#include <regex>
#include <string>
#include <typeinfo>

struct Unispg {
   protected:
        int refstart;
        int refend;
        int s;
        // int g_idx;

        int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
        int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
        //  Source and sink are also included.!!
        GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
   public:
        // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
        // b: all bundles on all strands: 0,1,2
   	    Unispg(int refstart_i=0, int refend_i=0, int s_i=0):refstart(refstart_i),refend(refend_i),s(s_i) { 
            no2gnode = new GPVec<CGraphnode>;
            no2gnode->setCapacity(8192);
            // fprintf(stderr, "no2gnode print test %d: ", no2gnode);
            CGraphnode* source=new CGraphnode(0,0,0);
            no2gnode -> Add(source);
            // fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(0)->nodeid);
            graphno = 1;
            edgeno = 0;
        }

        ~Unispg() {
            // no2gnode->Clear();
        }

        void Clear() {
            no2gnode -> Clear();
            no2gnode->setCapacity(8192);
        }


        void AddNode(int refstart_i, int refend_i, int start, int end, int node_id, float cov) {
            if (refstart_i == refstart && refend_i == refend) {
            CGraphnode* new_node=new CGraphnode(start, end, node_id, cov); // start,end,nodeno
            no2gnode -> Add(new_node);
            // fprintf(stderr, "Node size: %d\n ", no2gnode->Count());
            graphno += 1;
            }
        }

        void AddSink(int refstart_i, int refend_i) {
            if (refstart_i == refstart && refend_i == refend) {
            CGraphnode* sink=new CGraphnode(0,0,graphno);
            no2gnode -> Add(sink);
            // fprintf(stderr, "sink Node size: %d\n ", no2gnode->Count());
            // fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(graphno)->nodeid);
            graphno += 1;
            }
        }

        void AddEdge(int refstart_i, int refend_i, int head, int tail) {
            if (refstart_i == refstart && refend_i == refend) {
            // no2gnode[head];
            // fprintf(stderr, "no2gnode[0]->nodeid %d: ", no2gnode->Get(head));
            no2gnode->Get(head)->child.Add(tail);
            no2gnode->Get(tail)->parent.Add(head);

            // GBitVec childpat;
            // GBitVec parentpat;

            edgeno+=1;
            }
        }

        void PrintGraph() {
            { //DEBUG ONLY
            fprintf(stderr,"\tafter traverse:\n");
            for(int i=1;i<graphno-1;i++) {
                fprintf(stderr,"\tNode %d with parents:",i);
                for(int p=0;p<no2gnode->Get(i)->parent.Count();p++) fprintf(stderr," %d",no2gnode->Get(i)->parent[p]);
                fprintf(stderr," and children:");
                for(int c=0;c<no2gnode->Get(i)->child.Count();c++) fprintf(stderr," %d",no2gnode->Get(i)->child[c]);
            
                fprintf(stderr," %d(%d-%d)",i,no2gnode[0][i]->start,no2gnode[0][i]->end);
                fprintf(stderr,"\n");
            }
            }
        }

        int get_refstart() {
            return refstart;
        }
        int get_refend() {
            return refend;
        }
        int get_s() {
            return s;
        }
        int get_graphno() {
            return graphno;
        }
        int get_edgeno() {
            return edgeno;
        }
        GPVec<CGraphnode>* get_no2gnode() {
            return no2gnode;
        }

    // after the graph is created, I need to do 'predict transcripts for unstranded bundles here'
    /*****************************
     ** 'CPrediction': constructor
    **		predict transcripts for unstranded bundles here
    *****************************/

    // And then, I need to parse the graph here
    /*****************************
     ** 5. parse graph
    **    'process_refguides' & 'process_transfrags' & 'find_transcripts' & 'free_treepat'
    *****************************/
};

struct UnispgGp {
    protected:
        // //  Do not need
    // GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
        // //  Do not need
    // GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle : this might be the key for overalps
        // //  Do not need
        // int bno[2]={0,0};
        // int refstart;
        // int refend;
        // int gpSize[2];
        // GVec<int> graphnoGp[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
        // GVec<int> edgenoGp[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
        // int cgidx;
        GPVec<CGraphnode>* no2gnodeGp[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
        // GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
        // //  Do not need
        // GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
        // // GPVec<CMTransfrag> *mgt[2]; // merged super-transfrags

        // CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
        // GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
        // GVec<int> lastgpos[2];
    public:
        // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
        // b: all bundles on all strands: 0,1,2
        UnispgGp() { 
        // UnispgGp(int refstart_i, int refend_i): refstart(refstart_i), refend(refend_i) { 
            // cgidx = 0;
            for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
                int s=sno/2; // adjusted strand due to ignoring neutral strand
            //     gpSize[s] = 0;
            //     graphnoGp[s].setCapacity(8192);
            //     edgenoGp[s].setCapacity(8192);
            //     no2gnodeGp[s] = new GPVec<CGraphnode>[8192];
                no2gnodeGp[s] = new GPVec<CGraphnode>[100000];
            }
        }

        ~UnispgGp() {
            for(int i=0;i<2;i++) {
            // gpSize[i] = 0;
            // graphnoGp[i].Clear();
            // edgenoGp[i].Clear();
            delete [] no2gnodeGp[i];
            // no2gnodeGp[i] = new GPVec<CGraphnode>;
            // no2gnode[i] = NULL;
            };
        }

        void SetUnispgCapacity(int s) {
            fprintf(stderr, "**************************************\n");
            fprintf(stderr, "*********** SetUnispgCapacity ********\n");
            fprintf(stderr, "**************************************\n");
            // fprintf(stderr, "s: %d; capacity: %d\n", s, capacity);
            // gpSize[s] = 2;
            // graphnoGp[s].setCapacity(capacity);
            // edgenoGp[s].setCapacity(capacity);
            // graphnoGp[s].cAdd(0);
            // edgenoGp[s].cAdd(0);
            no2gnodeGp[s] = new GPVec<CGraphnode>[100000];
        }

        // void SetRefStartEnd(int refstart_i, int refend_i) {
        //     refstart = refstart_i;
        //     refend = refend_i;
        // }

        // void UpdateRefStartEnd(int refstart_i, int refend_i) {
        //     refstart = refstart_i;
        //     refend = refend_i;
        // }

        void AddGraph(int fidx, int s, int cgidx, GPVec<CGraphnode>* no2gnode) {
                // GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
                // GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
                // GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i


            // fprintf(stderr, "* uni_splice_graph.refstart: %d \n", uni_splice_graph -> get_refstart());
            // fprintf(stderr, "* uni_splice_graph.refend: %d \n", uni_splice_graph -> get_refend());
            // fprintf(stderr, "* uni_splice_graph.s: %d \n", uni_splice_graph -> get_s());
            // fprintf(stderr, "* uni_splice_graph.g_idx: %d \n", uni_splice_graph -> get_g_idx());
            // fprintf(stderr, "* uni_splice_graph.get_graphno: %d \n", uni_splice_graph -> get_graphno());
            // fprintf(stderr, "* uni_splice_graph.get_edgeno: %d \n", uni_splice_graph -> get_edgeno());

            // int refstart;
            // int refend;
            // int s;
            // int g_idx;

            // int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
            // int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
            // GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
            if (fidx == 0) {
                fprintf(stderr, "*****************************\n");
                fprintf(stderr, "*********** AddGraph ********\n");
                fprintf(stderr, "*****************************\n");
                // This is the first unispg. Simply add copy it into UnispgGp
                // graphnoGp[s].Add(graphno);
                // edgenoGp[s].Add(edgeno);
                // gpSize[s] = no2gnode->Count();
                // no2gnodeGp[s][b].Add(no2gnode->Get(0));

                for (int i=0; i<no2gnode->Count(); i++) {
                    no2gnodeGp[s][cgidx].Add(no2gnode->Get(i));
                    // cgidx = cgidx+1;
                    fprintf(stderr, "Inside!! cgidx: %d \n", cgidx);
                }   
            } else {
                fprintf(stderr, "************************************\n");
                fprintf(stderr, "*********** Add more graph! ********\n");
                fprintf(stderr, "************************************\n");
                // // Need to check overlapping & new graph creation
                // for(int sno=0; sno<3; sno+=2) { // skip neutral bundles -> those shouldn't have junctions
                //     int s=sno/2; // adjusted strand due to ignoring neutral strand
                //     // fprintf(stderr, "2 bundle[sno].Count(): %d\n", bundle[sno].Count());
                //     unispg_gp->SetUnispgCapacity(s, bundle[sno].Count());
                //     for(int b=0;b<bundle[sno].Count();b++) {
                //         // Unispg unispg = new Unispg();
                //         unispg_gp->AddGraph(fidx, s, b, graphno[s][b], edgeno[s][b], no2gnode[s]+b);
                //     }
                // }
            }

            // int g_idx = 0;
            // for (int i=0; i<no2gnode->Count(); i++) {
            //     no2gnodeGp[s][g_idx].Add(no2gnode->Get(i));
            //     g_idx++;
            // }


            // gpSize[uni_splice_graph->get_s()] = uni_splice_graph->get_g_idx()+1;
            // fprintf(stderr, "&&& uni_splice_graph->get_g_idx()+1: %d\n", uni_splice_graph->get_g_idx()+1);
            // fprintf(stderr, "&&& gpSize[uni_splice_graph->get_s()]: %d\n", gpSize[uni_splice_graph->get_s()]);
            // We need to reset graph_idx when (1)the new strands     
        }


        void Clear() {
            // fprintf(stderr, "**** Start Clearing !!!! \n ");
            for(int i=0;i<2;i++) {
            // gpSize[i] = 0;
            // graphnoGp[i].Clear();
            // graphnoGp[i].setCapacity(8192);
            // edgenoGp[i].Clear();
            // edgenoGp[i].setCapacity(8192);
            delete [] no2gnodeGp[i];
            no2gnodeGp[i] = new GPVec<CGraphnode>[8192];
            // no2gnodeGp[i]->setCapacity(8192);
            };
        }

// Need to be modifed!!!
        void PrintGraphGp() {
            fprintf(stderr, "*********************************\n");
            fprintf(stderr, "*********** PrintGraphGp ********\n");
            fprintf(stderr, "*********************************\n");

            { // DEBUG ONLY
                printTime(stderr);
                // for(int s=0;s<2;s++) {
                //     fprintf(stderr, "\n\tThere are %d stranded[%d] graphs\n", gpSize[s],int(2*s));
                //     // if () {

                //     // }
                //     for(int b=0;b<gpSize[s];b++) {
                //         fprintf(stderr, ">>>>>>> 1-2 gpSize[%d]: %d\n", s, gpSize[s]);
                //         fprintf(stderr, ">>>>>>> 1-2 graphnoGp[%d][%d]: %d\n", s, b, graphnoGp[s][b]);
                //         if(graphnoGp[s][b]) {
                //             GStr pat;
                //             fprintf(stderr,"\t\tGraph[%d][%d] with %d nodes and %d edges :",int(2*s),b,graphnoGp[s][b],edgenoGp[s][b]);
                //             for(int nd=1;nd<graphnoGp[s][b]-1;nd++) {
                //                 fprintf(stderr, "&&& nd: %d\n", nd);
                //                 fprintf(stderr, "&&& no2gnodeGp[s][b][nd]->start: %d\n", no2gnodeGp[s][b][nd]->start);
                //                 fprintf(stderr, "&&& no2gnodeGp[s][b][nd]->end: %d\n", no2gnodeGp[s][b][nd]->end);
                //                 fprintf(stderr," %d(%d-%d)",nd,no2gnodeGp[s][b][nd]->start,no2gnodeGp[s][b][nd]->end);
                //             }
                //         } 
                //         fprintf(stderr,"\n");
                //     }
                // }

            //     for(int s=0;s<2;s++) {
    		// 	fprintf(stderr, "There are %d stranded[%d] graphs\n",bno[s],int(2*s));
			// 	// fprintf(stderr, "1 bundle[sno].Count(): %d\n", bno[s]);
    		// 	for(int b=0;b<bno[s];b++) {
    		// 		if(graphno[s][b]) {
    		// 			GStr pat;
    		// 			fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges with lastgpos=%d:",int(2*s),b,graphno[s][b],edgeno[s][b],lastgpos[s][b]);
    		// 			for(int nd=1;nd<graphno[s][b]-1;nd++)
    		// 				fprintf(stderr," %d(%d-%d)",nd,no2gnode[s][b][nd]->start,no2gnode[s][b][nd]->end);
    		// 			fprintf(stderr,"\n");
    		// 			print_pattern(tr2no[s][b],pat,graphno[s][b]);
    		// 		}
    		// 	}
    		// }
            }
        }

        // int get_refstart () {
        //     return refstart;
        // }
        // int get_refend () {
        //     return refend;
        // }
        // GVec<int>* get_graphnoGp () {
        //     return graphnoGp;
        // }
        // GVec<int>* get_edgenoGp () {
        //     return edgenoGp;
        // }
        GPVec<CGraphnode>** get_no2gnodeGp () {
            return no2gnodeGp;
        }
        // int* get_gpSize() {
        //     return gpSize;
        // }
};


#endif