#ifndef __MULTIST_H__
#define __MULTIST_H__
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

struct UniSpliceGraph {
   protected:
      int refstart;
	   int refend;
      int s;
      int g_idx;

      int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
      int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
      GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
   public:
   // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
   // b: all bundles on all strands: 0,1,2
   	UniSpliceGraph(int refstart_i=0, int refend_i=0, int s_is=0, int g_idx_i=0):refstart(refstart_i),refend(refend_i),s(s_is),g_idx(g_idx_i) { 
    		no2gnode = new GPVec<CGraphnode>;
         // fprintf(stderr, "no2gnode print test %d: ", no2gnode);
         CGraphnode* source=new CGraphnode(0,0,0);
         no2gnode -> Add(source);
         fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(0)->nodeid);
         graphno = 0;
         edgeno = 0;
      }

      void AddNode(int refstart_i, int refend_i, int start, int end, int node_id) {
         if (refstart_i == refstart && refend_i == refend) {
            CGraphnode* new_node=new CGraphnode(start, end, node_id); // start,end,nodeno
            no2gnode -> Add(new_node);
            fprintf(stderr, "Node size: %d\n ", no2gnode->Count());
            graphno += 1;
         }
      }

      void AddSink(int refstart_i, int refend_i) {
         if (refstart_i == refstart && refend_i == refend) {
            CGraphnode* sink=new CGraphnode(0,0,graphno+1);
            no2gnode->Add(sink);
            fprintf(stderr, "sink Node size: %d\n ", no2gnode->Count());
            fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(graphno+1)->nodeid);
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
            fprintf(stderr,"after traverse:\n");
            for(int i=1;i<graphno+1;i++) {
               fprintf(stderr,"Node %d with parents:",i);
               for(int p=0;p<no2gnode->Get(i)->parent.Count();p++) fprintf(stderr," %d",no2gnode->Get(i)->parent[p]);
               fprintf(stderr," and children:");
               for(int c=0;c<no2gnode->Get(i)->child.Count();c++) fprintf(stderr," %d",no2gnode->Get(i)->child[c]);
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
      int get_g_idx() {
         return g_idx;
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

struct UniSpliceGraphGp {
   protected:
      // //  Do not need
   	// GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
      // //  Do not need
   	// GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle : this might be the key for overalps
      // //  Do not need
    	// int bno[2]={0,0};
      GVec<int> graphnoGp[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
      GVec<int> edgenoGp[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
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
   	UniSpliceGraphGp() { 
         for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
            int s=sno/2; // adjusted strand due to ignoring neutral strand
            graphnoGp[s].setCapacity(4096);
            edgenoGp[s].setCapacity(4096);

            no2gnodeGp[s] = new GPVec<CGraphnode>;
         }
      }

      void AddGraph(UniSpliceGraph* uni_splice_graph) {
         fprintf(stderr, "* uni_splice_graph.refstart: %d \n", uni_splice_graph -> get_refstart());
         fprintf(stderr, "* uni_splice_graph.refend: %d \n", uni_splice_graph -> get_refend());
         fprintf(stderr, "* uni_splice_graph.s: %d \n", uni_splice_graph -> get_s());
         fprintf(stderr, "* uni_splice_graph.g_idx: %d \n", uni_splice_graph -> get_g_idx());
         fprintf(stderr, "* uni_splice_graph.get_graphno: %d \n", uni_splice_graph -> get_graphno());
         fprintf(stderr, "* uni_splice_graph.get_edgeno: %d \n", uni_splice_graph -> get_edgeno());

      //    int refstart;
	   // int refend;
      // int s;
      // int g_idx;

      // int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
      // int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
      // GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
         graphnoGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()] = uni_splice_graph->get_graphno();

         // fprintf(stderr, "* After graphno[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()]: \n");
         
         edgenoGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()] = uni_splice_graph->get_edgeno();
         // fprintf(stderr, "* After graphno[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()]: \n");

         // fprintf(stderr, "* graphno[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()]: %d \n", graphnoGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()]);
         // fprintf(stderr, "* edgeno[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()]: %d \n", edgenoGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()]);
         // no2gnodeGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()] = no2gnode;

         // CGraphnode* source=new CGraphnode(0,0,0);
         no2gnodeGp[uni_splice_graph->get_s()] = uni_splice_graph->get_no2gnode();
         // We need to reset graph_idx when (1)the new strands   
      }

      void Clear() {
         for(int i=0;i<2;i++) {
            graphnoGp[i].Clear();
            edgenoGp[i].Clear();
            graphnoGp[i].setCapacity(1024);
            edgenoGp[i].setCapacity(1024);
            // delete no2gnode[i];
            no2gnodeGp[i]->Clear();
            no2gnodeGp[i] = new GPVec<CGraphnode>;
            // no2gnode[i] = NULL;
         };
         // fprintf(stderr, "Done cleaning!! \n");
	   }

      void PrintGraphGp() {
         fprintf(stderr, "*********************************\n");
         fprintf(stderr, "*********** PrintGraphGp ********\n");
         fprintf(stderr, "*********************************\n");

         { // DEBUG ONLY
            printTime(stderr);
            for(int s=0;s<2;s++) {
               fprintf(stderr, "There are %d stranded[%d] graphs\n",no2gnodeGp[s]->Count(),int(2*s));
               for(int b=0;b<no2gnodeGp[s]->Count();b++) {
                  if(graphnoGp[s][b]) {
                     GStr pat;
                     fprintf(stderr,"Graph[%d][%d] with %d nodes and %d edges :",int(2*s),b,graphnoGp[s][b],edgenoGp[s][b]);
                     for(int nd=1;nd<graphnoGp[s][b]-1;nd++)
                        fprintf(stderr," %d(%d-%d)",nd,no2gnodeGp[s][b][nd]->start,no2gnodeGp[s][b][nd]->end);
                     fprintf(stderr,"\n");
                  }
               }
            }
         }
      }
      GVec<int>* get_graphnoGp () {
         return graphnoGp;
      }
      GVec<int>* get_edgenoGp () {
         return edgenoGp;
      }
      GPVec<CGraphnode>** get_no2gnodeGp () {
         return no2gnodeGp;
      }
};



      

class DOTReader;
class DOTWriter;

// class DOTRecord {
//    friend class DOTReader;
//    friend class DOTWriter;
   
//    public:
//       UniSpliceGraph* uni_Splice_Graph;
//       DOTRecord() {};

//       void init() {
//          clear();
//       }

//       void clear() {
//       }

//       ~DOTRecord() {
//          clear();
//       }

//       void parse_error(const char* s) {
//          GError("DOT parsing error: %s\n", s);
//       }
// };

class DOTReader {
   public:
      const char* fname;
      ifstream ifile_dot;
      bool dotopen(const char* filename) {
	      if (ifile_dot.is_open()) {
            return ifile_dot.is_open();
         } else {
            GError("Error: could not open the dot file %s \n",filename);
         }
      }

      DOTReader(const char* fn) {
         ifile_dot.open(fn);
      }

      const char* fileName() {
         return fname;
      }

      void dotclose() {
         ifile_dot.close();
         // if (hts_file) {
         //       if (hdr!=NULL) sam_hdr_destroy(hdr);
         //       hdr=NULL;
         //    hts_close(hts_file);
         //    hts_file=NULL;
         // }
      }

      ~DOTReader() {
         dotclose();
         GFREE(fname);
      }
      
      //the caller has to FREE the created UniSpliceGraph
      UniSpliceGraph* next() {
         // Three target Vector to be created.
         // GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
         // GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
         // GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i



         // Example from: https://www.codegrepper.com/code-examples/cpp/remove+all+spaces+from+string+c%2B%2B
         string line;  
         if (getline(ifile_dot, line)) {
            // line.erase(remove(line.begin(), line.end(), ' '), line.end());
            fprintf(stderr, "Before parsing line: %s \n", line.c_str());
            // if (regex_match(line, regex("(strict\\s+)(digraph\\s+)([0-9]*)(_)([0-9]*)(_)([0-9]*\\s+)(.*)(->)(.*)(\\[label=)(.*)(\\];)"))) {
            // regex rgx("(\\w+)->(\\w+)\\[label=(\\w+)\\];");

            string delimiter = "{";
            size_t pos = 0;
            string token;



            // for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
            //    int s=sno/2; // adjusted strand due to ignoring neutral strand
            //    no2gnode[s] = new GPVec<CGraphnode>;
            // }


            int refstart = 0;
            int refend = 0;
            int s = 0;
            int g_idx = 0;
            bool add_sink = true;
            while ((pos = line.find(delimiter)) != string::npos) {
               token = line.substr(0, pos);
               // std::cout << token << std::endl;
               // fprintf(stderr, "token: %s\n", token.c_str());
               regex title_rgx("strict\\s+digraph\\s+(\\w+)_(\\w+)_(\\w+)_(\\w+)");
               smatch match;
               // fprintf(stderr, "line: %s \n", line.c_str());
               if (regex_search(token, match, title_rgx)) {
                  refstart = stoi(match[1]);
                  refend = stoi(match[2]);
                  s = stoi(match[3]);
                  g_idx = stoi(match[4]);
                  fprintf(stderr, "ref : %d - %d, %d, %d\n", refstart, refend, s, g_idx);
                  // fprintf(stderr, "refstart: %d\n", refend);
                  // fprintf(stderr, "s: %d\n", s);
                  // fprintf(stderr, "g_idx: %d\n", g_idx);
               }

               line.erase(0, pos + delimiter.length());
               break;
            }

            UniSpliceGraph* uni_splice_graph = new UniSpliceGraph(refstart, refend, s, g_idx);

            delimiter = ";";
            pos = 0;
            token = "";

            while ((pos = line.find(delimiter)) != string::npos) {
               int node_id = 0;
               int start = 0;
               int end = 0;
               int cov = 0;
               int head = 0;
               int tail = 0;
               token = line.substr(0, pos);
               // std::cout << token << std::endl;
               // fprintf(stderr, "token: %s\n", token.c_str());

               regex node_rgx("(\\w+)\\[+start=(\\w+)\\s+end=(\\w+)\\s+cov=(\\w+)");
               regex edge_rgx("(\\w+)->(\\w+)");
               smatch match;
               if (regex_search(token, match, node_rgx)) {
                  node_id = stoi(match[1]);
                  start = stoi(match[2]);
                  end = stoi(match[3]);
                  cov = stoi(match[4]);
                  // fprintf(stderr, "node : %d\n", node);
                  // fprintf(stderr, "start : %d\n", start);
                  // fprintf(stderr, "end : %d\n", end);
                  // fprintf(stderr, "cov : %d\n", cov);
                  uni_splice_graph->AddNode(refstart, refend, start, end, node_id);
               }
               if (regex_search(token, match, edge_rgx)) {
                  // add sink
                  if (add_sink) {
                     uni_splice_graph->AddSink(refstart, refend);
                     add_sink = false;
                  }
                  head = stoi(match[1]);
                  tail = stoi(match[2]);
                  // fprintf(stderr, "head : %d\n", head);
                  // fprintf(stderr, "tail : %d\n", tail);
                  uni_splice_graph->AddEdge(refstart, refend, head, tail);
               }
               line.erase(0, pos + delimiter.length());

               // no2gnode
            }
            uni_splice_graph -> PrintGraph();

            return uni_splice_graph;
         } else {    
            return NULL;        
         }   
      }

      bool next(UniSpliceGraph& rec) {
         return false;
      }
};


//basic BAM/SAM/CRAM writer class
// limitations: cannot add new reference sequences, just new alignments to
//  existing reference sequences;
class DOTWriter {
 public:
   void create() {
     create();
   }

   DOTWriter() {
   }

   ~DOTWriter() {
   }
};




/*********************************
 * bundle-related data structure.
 *********************************/
// // bundle data structure, holds all data needed for
// // infering transcripts from a bundle
// enum BundleStatus {
// 	BUNDLE_STATUS_CLEAR=0, //available for loading/prepping
// 	BUNDLE_STATUS_LOADING, //being prepared by the main thread (there can be only one)
// 	BUNDLE_STATUS_READY //ready to be processed, or being processed
// };

// struct CBundle {
// 	int len;
// 	float cov;
// 	float multi;
// 	int startnode;  // id of start node in bundle of same strand
// 	int lastnodeid; // id of last node added to bundle
// 	CBundle(int _len=0, float _cov=0, float _multi=0, int _start=-1, int _last=-1):
// 			len(_len),cov(_cov),multi(_multi), startnode(_start),lastnodeid(_last) {}
// };

// struct CBundlenode:public GSeg {
// 	float cov;
// 	int bid; // bundle node id in bnode -> to easy retrieve it
// 	CBundlenode *nextnode; // next node in the same bundle
// 	CBundlenode(int rstart=0, int rend=0, float _cov=0, int _bid=-1, CBundlenode *_nextnode=NULL):GSeg(rstart, rend),
// 			cov(_cov),bid(_bid),nextnode(_nextnode) {}
// };

// { // DEBUG ONLY
//     fprintf(stderr,"process refguides for s=%d b=%d edgeno=%d gno=%d lastgpos=%d guidescount=%d\n",s,b,edgeno[s][b],graphno[s][b],lastgpos[s][b],guides.Count());
//     fprintf(stderr,"There are %d nodes for graph[%d][%d]:\n",graphno[s][b],s,b);
//     for(int i=0;i<graphno[s][b];i++) {
//         fprintf(stderr,"%d (%d-%d): %f len=%d cov=%f",i,no2gnode[s][b][i]->start,no2gnode[s][b][i]->end,no2gnode[s][b][i]->cov,no2gnode[s][b][i]->len(),no2gnode[s][b][i]->cov/no2gnode[s][b][i]->len());
//         fprintf(stderr," parents:");
//         for(int j=0;j<no2gnode[s][b][i]->parent.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->parent[j]);
//         fprintf(stderr," trf=");
//         for(int j=0;j<no2gnode[s][b][i]->trf.Count();j++) fprintf(stderr," %d",no2gnode[s][b][i]->trf[j]);
//         fprintf(stderr,"\n");
//     }
// }


// edgeno
// graphno
// lastgpos
// guides
// no2gnode

/*********************************
 * graph-related data structure.
 *********************************/

// struct CGraphinfo {
// 	int ngraph;
// 	int nodeno;
// 	CGraphinfo(int ng=-1,int nnode=-1):ngraph(ng),nodeno(nnode){}
// };

// struct CGraphnode:public GSeg {
// 	int nodeid;
// 	float cov;
// 	float capacity; // sum of all transcripts abundances exiting and through node
// 	float rate; // conversion rate between in and out transfrags of node
// 	//float frag; // number of fragments included in node
// 	GVec<int> child;
// 	GVec<int> parent;
// 	GBitVec childpat;
// 	GBitVec parentpat;
// 	GVec<int> trf; // transfrags that pass the node
// 	bool hardstart:1; // verified/strong start
// 	bool hardend:1;	// verified/strong end
// 	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
// 	CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0):GSeg(s,e),
// 			nodeid(id),cov(nodecov),capacity(cap),rate(r),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){}
// };

// struct DOTInputRecord {
//    UniSpliceGraph* brec;
//    int fidx; //index in files and readers

// 	DOTInputRecord(UniSpliceGraph* b=NULL, int i=0):brec(b),fidx(i) {}
// 	~DOTInputRecord() {
// 		delete brec;
// 	}
// };


struct DOTInputFile {
   protected:
      UniSpliceGraph* rec;
      // UniSpliceGraphGp* uni_splice_graphs;
   public:
      DOTReader* reader;
      GStr file; //same order
      GStr tmpfile; //all the temp files created by this
      DOTInputFile():rec(NULL), reader(), file(), tmpfile() {
         fprintf(stderr, "&& DOTInputFile Initialization  :\n");
         // uni_splice_graphs = new UniSpliceGraphGp();
         // fprintf(stderr, "Initialization uni_splice_graphs.graph_idx %d :\n", uni_splice_graphs->get_graphno());
      }
      // void Add(const char* fn);
      int count() { return 1; }
      bool start(const char* fn); //open all files, load 1 record from each
      UniSpliceGraph* next();
      void stop();
};


#endif