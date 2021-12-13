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

struct UniSpliceGraph {
   protected:
      //  Do not need
   	GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
      //  Do not need
   	GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle : this might be the key for overalps
      //  Do not need
    	int bno[2]={0,0};


      int refstart;
	   int refend;
      int s;
      int g_idx;
      GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
      GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
      GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i




      GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
      //  Do not need
      GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    	// GPVec<CMTransfrag> *mgt[2]; // merged super-transfrags

      CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
      GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
      GVec<int> lastgpos[2];
   public:
   // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
   // b: all bundles on all strands: 0,1,2
   	UniSpliceGraph():refstart(0), refend(0), s(NULL), g_idx(NULL) { 
         //    // build graph structure
         //    for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions

         //       int s=sno/2; // adjusted strand due to ignoring neutral strand
         //       //char strnd='-';
         //       //if(s) strnd='+';

         //       bundle2graph[s]=NULL;
         //       if(bnode[sno].Count()) bundle2graph[s]=new GVec<CGraphinfo>[bnode[sno].Count()];
         //       transfrag[s]=NULL;
         //       // mgt[s]=NULL;
         //       no2gnode[s]=NULL;
         //       tr2no[s]=NULL;
         //       gpos[s]=NULL;

         //       if(bundle[sno].Count()) {
         //          transfrag[s]=new GPVec<CTransfrag>[bundle[sno].Count()]; // for each bundle I have a graph ? only if I don't ignore the short bundles
         //          // mgt[s]=new GPVec<CMTransfrag>[bundle[sno].Count()];
         //          no2gnode[s]=new GPVec<CGraphnode>[bundle[sno].Count()];

         //          gpos[s]=new GIntHash<int>[bundle[sno].Count()];


         //          GCALLOC(tr2no[s],bundle[sno].Count()*sizeof(CTreePat *));
         //          bno[s]=bundle[sno].Count();

         //          for(int b=0;b<bundle[sno].Count();b++) {
         //             graphno[s].cAdd(0);
         //             edgeno[s].cAdd(0);
         //             lastgpos[s].cAdd(0);

         //             /*
         //             { // DEBUG ONLY
         //             if(bundle[sno][b]->nread) {
         //                fprintf(stderr,"proc bundle[%d][%d] %f/%f is %f len=%d\n",sno,b,bundle[sno][b]->multi,bundle[sno][b]->nread,(float)bundle[sno][b]->multi/bundle[sno][b]->nread,bundle[sno][b]->len);
         //             } }
         //             */

         //             // after the graph is created, I need to do 'predict transcripts for unstranded bundles here'
         //             /*****************************
         //              ** I need to come up with my own `create_unisplicegraph`
         //              *****************************/
         //             // create graph
         //             GArray<GEdge> unused;
         //             graphno[s][b]=create_graph(refstart,s,b,bundle[sno][b],bnode[sno],junction,ejunction,
         //                   bundle2graph,no2gnode,transfrag,gpos,NULL,edgeno[s][b],lastgpos[s][b],unused,refend);

         //             if(graphno[s][b]) tr2no[s][b]=construct_treepat(graphno[s][b],gpos[s][b],transfrag[s][b]);
         //             else tr2no[s][b]=NULL;

         //             // for(int i=0; i<transfrag[s][b].Count();i++) {
         //             //    CMTransfrag *tr=new CMTransfrag(transfrag[s][b][i]); // transfrags that are created in the graph process can not be associated with any read
         //             //    mgt[s][b].Add(tr);
         //             // }

         //          }
         //       }
         //    }
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
      

class DOTReader;
class DOTWriter;

class DOTRecord {
   friend class DOTReader;
   friend class DOTWriter;
   
   public:
      UniSpliceGraph uni_Splice_Graph;
      DOTRecord();

      void init() {
         clear();
      }

      void clear() {
      }

      ~DOTRecord() {
         clear();
      }

      void parse_error(const char* s) {
         GError("DOT parsing error: %s\n", s);
      }
};

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
      
      //the caller has to FREE the created DOTRecord
      DOTRecord* next() {

         // Three target Vector to be created.
         GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
         GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
         GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i



         // Example from: https://www.codegrepper.com/code-examples/cpp/remove+all+spaces+from+string+c%2B%2B
         string line;
         getline(ifile_dot, line);         
         // line.erase(remove(line.begin(), line.end(), ' '), line.end());
         fprintf(stderr, "line: %s \n", line.c_str());
         // if (regex_match(line, regex("(strict\\s+)(digraph\\s+)([0-9]*)(_)([0-9]*)(_)([0-9]*\\s+)(.*)(->)(.*)(\\[label=)(.*)(\\];)"))) {
         // regex rgx("(\\w+)->(\\w+)\\[label=(\\w+)\\];");

         string delimiter = "{";
         size_t pos = 0;
         string token;



         for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
            int s=sno/2; // adjusted strand due to ignoring neutral strand
    		   no2gnode[s]=NULL;
         }


         int refstart = 0;
         int refend = 0;
         int s = 0;
         int g_idx = 0;
         while ((pos = line.find(delimiter)) != string::npos) {
            token = line.substr(0, pos);
            // std::cout << token << std::endl;
            fprintf(stderr, "token: %s\n", token.c_str());
            regex title_rgx("strict\\s+digraph\\s+(\\w+)_(\\w+)_(\\w+)_(\\w+)");
            smatch match;
            fprintf(stderr, "line: %s \n", line.c_str());
            if (regex_search(token, match, title_rgx)) {
               refstart = stoi(match[1]);
               refend = stoi(match[2]);
               s = stoi(match[3]);
               g_idx = stoi(match[4]);
               fprintf(stderr, "refstart: %d\n", refstart);
               fprintf(stderr, "refstart: %d\n", refend);
               fprintf(stderr, "s: %d\n", s);
               fprintf(stderr, "g_idx: %d\n", g_idx);
            }

            line.erase(0, pos + delimiter.length());
            break;
         }


    		// no2gnode[s]=new GPVec<CGraphnode>[bundle[sno].Count()];
         delimiter = ";";
         pos = 0;
         token = "";

         // Add the source first.
         // CGraphnode* source=new CGraphnode(0,0,0);
         // no2gnode[s][g_idx].Add(source);
	      // CGraphnode* sink=new CGraphnode();



         while ((pos = line.find(delimiter)) != string::npos) {
            int node = 0;
            int start = 0;
            int end = 0;
            int cov = 0;
            int head = 0;
            int tail = 0;
            token = line.substr(0, pos);
            // std::cout << token << std::endl;
            fprintf(stderr, "token: %s\n", token.c_str());

            regex node_rgx("(\\w+)\\[+start=(\\w+)\\s+end=(\\w+)\\s+cov=(\\w+)");
            regex edge_rgx("(\\w+)->(\\w+)");
            smatch match;
            if (regex_search(token, match, node_rgx)) {
               node = stoi(match[1]);
               start = stoi(match[2]);
               end = stoi(match[3]);
               cov = stoi(match[4]);
               fprintf(stderr, "node : %d\n", node);
               fprintf(stderr, "start : %d\n", start);
               fprintf(stderr, "end : %d\n", end);
               fprintf(stderr, "cov : %d\n", cov);
            }
            if (regex_search(token, match, edge_rgx)) {
               head = stoi(match[1]);
               tail = stoi(match[2]);
               fprintf(stderr, "head : %d\n", head);
               fprintf(stderr, "tail : %d\n", tail);
            }
            line.erase(0, pos + delimiter.length());

            // no2gnode
         }

         return NULL;
      }
      bool next(DOTRecord& rec) {
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
//    DOTRecord* brec;
//    int fidx; //index in files and readers

// 	DOTInputRecord(DOTRecord* b=NULL, int i=0):brec(b),fidx(i) {}
// 	~DOTInputRecord() {
// 		delete brec;
// 	}
// };


struct DOTInputFile {
   protected:
      DOTRecord* rec;
   public:
      DOTReader* reader;
      GStr file; //same order
      GStr tmpfile; //all the temp files created by this
      // GList<DOTInputRecord> recs; //next record for each
      DOTInputFile():rec(NULL), reader(), file(), tmpfile() { }
      // void Add(const char* fn);
      int count() { return 1; }
      bool start(const char* fn); //open all files, load 1 record from each
      DOTRecord* next();
      void stop();
};


#endif
