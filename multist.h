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
    GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
    GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
    GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
    // GVec<int> trnumber[2]; // how many transfrags are on a strand s, in a graph g -> I can find out this from transfrag[s][g].Count()
    // int ngraph[2]={0,0};   // how many graphs are in each strand: negative (0), or positive(1) -> keep one for each bundle
    GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
    GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
    CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
    GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
    GVec<int> lastgpos[2];
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
         while ((pos = line.find(delimiter)) != string::npos) {
            token = line.substr(0, pos);
            // std::cout << token << std::endl;
            fprintf(stderr, "token: %s\n", token.c_str());
            regex title_rgx("strict\\s+digraph\\s+(\\w+)_(\\w+)_(\\w+)");
            smatch match;
            fprintf(stderr, "line: %s \n", line.c_str());
            if (regex_search(token, match, title_rgx)) {
               string refstart = match[1];
               string s = match[2];
               string b = match[2];
               fprintf(stderr, "refstart: %s\n", refstart.c_str());
               fprintf(stderr, "s: %s\n", s.c_str());
               fprintf(stderr, "b: %s\n", b.c_str());
            }

            line.erase(0, pos + delimiter.length());
            break;
         }

         delimiter = ";";
         pos = 0;
         token = "";
         while ((pos = line.find(delimiter)) != string::npos) {
            token = line.substr(0, pos);
            // std::cout << token << std::endl;
            fprintf(stderr, "token: %s\n", token.c_str());

            regex node_rgx("(\\w+)\\[+start=(\\w+)\\s+end=(\\w+)\\s+cov=(\\w+)");
            regex edge_rgx("(\\w+)->(\\w+)");
            smatch match;
            if (regex_search(token, match, node_rgx)) {
               string node = match[1];
               string start = match[2];
               string end = match[3];
               string cov = match[4];
               fprintf(stderr, "node : %s\n", node.c_str());
               fprintf(stderr, "start : %s\n", start.c_str());
               fprintf(stderr, "end : %s\n", end.c_str());
               fprintf(stderr, "cov : %s\n", cov.c_str());
            }
            if (regex_search(token, match, edge_rgx)) {
               string head = match[1];
               string tail = match[2];
               fprintf(stderr, "head : %s\n", head.c_str());
               fprintf(stderr, "tail : %s\n", tail.c_str());
            }

            line.erase(0, pos + delimiter.length());
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
