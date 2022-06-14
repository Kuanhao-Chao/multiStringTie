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
#include "unispg.h"
#include <fstream>
#include <regex>
#include <string>
#include <typeinfo>

// struct CGraphnodeUnispg:public GSeg {
//     int sample_num = 0;
// 	int old_graph_id = -1;
// 	int old_node_id = -1;
// 	int nodeid;
// 	// Samples having this node
//     GVec<bool>* is_passed_s;
// 	// Node coverage for each sample
//     GVec<float>* cov_s;
// 	// Node capacity for each sample
//     GVec<float>* capacity_s;

// 	GVec<int> child;
// 	GBitVec childpat;
// 	GBitVec parentpat;
// 	GVec<int> trf; // transfrags that pass the node
// 	bool hardstart:1; // verified/strong start
// 	bool hardend:1;	// verified/strong end
// 	//CGraphnode(int s=0,int e=0,unsigned int id=MAX_NODE,float nodecov=0,float cap=0,float r=0,float f=0):GSeg(s,e),nodeid(id),cov(nodecov),capacity(cap),rate(r),frag(f),child(),parent(),childpat(),parentpat(),trf(){}
// 	// CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0, int old_graph_id_i=0, int old_node_id_i=0, unsigned int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, GVec<float>* cov_s_i=NULL, GVec<float>* capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0):GSeg(s,e),sample_num(sample_num_i), old_graph_id(old_graph_id_i), old_node_id(old_node_id_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i),child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){
// 	CGraphnodeUnispg(int sample_num_i=0, int s=0,int e=0, int id=MAX_NODE, GVec<bool>* is_passed_s_i=NULL, GVec<float>* cov_s_i=NULL, GVec<float>* capacity_s_i=NULL, bool is_passed=false, float cov=0, float capacity=0,float r=0, bool set_g_n_idx=false, int g_idx=-1, int n_idx=-1):GSeg(s,e),sample_num(sample_num_i), nodeid(id),is_passed_s(is_passed_s_i),cov_s(cov_s_i),capacity_s(capacity_s_i), child(),parent(),childpat(),parentpat(),trf(),hardstart(false),hardend(false){

// 		fprintf(stderr, "		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
// 		fprintf(stderr, "		^^^ Creating graphnode id: %d (%u - %u) \n", id, s, e);
// 		fprintf(stderr, "		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
// 		fprintf(stderr, "		^^ is_passed: %d\n", is_passed);
// 		fprintf(stderr, "		^^ cov      : %f\n", cov);
// 		fprintf(stderr, "		^^ capacity : %f\n", capacity);
// 		is_passed_s->cAdd(is_passed);
// 		cov_s->cAdd(cov);
// 		capacity_s->cAdd(capacity);
// 		if (set_g_n_idx) {
// 			old_graph_id = g_idx;
// 			old_node_id = n_idx;
// 		}
//     }

// 	void set_nodeid(int new_node_id) {
// 		nodeid = new_node_id;
// 	}
// };


// struct UniSpliceGraph {
//    protected:
//       int refstart;
// 	   int refend;
//       int s;
//       int g_idx;

//       int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
//       int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
//       //  Source and sink are also included.!!
//       GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
//    public:
//       // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
//       // b: all bundles on all strands: 0,1,2
//    	UniSpliceGraph(int refstart_i=0, int refend_i=0, int s_is=0, int g_idx_i=0):refstart(refstart_i),refend(refend_i),s(s_is),g_idx(g_idx_i) { 
//     		no2gnode = new GPVec<CGraphnode>;
//          no2gnode->setCapacity(8192);
//          // fprintf(stderr, "no2gnode print test %d: ", no2gnode);
//          CGraphnode* source=new CGraphnode(0,0,0);
//          no2gnode -> Add(source);
//          // fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(0)->nodeid);
//          graphno = 1;
//          edgeno = 0;
//       }

//       ~UniSpliceGraph() {
//          // no2gnode->Clear();
//       }

//       void Clear() {
//          no2gnode -> Clear();
//          no2gnode->setCapacity(8192);
//       }

   
//       void AddNode(int refstart_i, int refend_i, int start, int end, int node_id, float cov) {
//          if (refstart_i == refstart && refend_i == refend) {
//             CGraphnode* new_node=new CGraphnode(start, end, node_id, cov); // start,end,nodeno
//             no2gnode -> Add(new_node);
//             // fprintf(stderr, "Node size: %d\n ", no2gnode->Count());
//             graphno += 1;
//          }
//       }

//       void AddSink(int refstart_i, int refend_i) {
//          if (refstart_i == refstart && refend_i == refend) {
//             CGraphnode* sink=new CGraphnode(0,0,graphno);
//             no2gnode -> Add(sink);
//             // fprintf(stderr, "sink Node size: %d\n ", no2gnode->Count());
//             // fprintf(stderr, "no2gnode[0]->nodeid %d: \n", no2gnode->Get(graphno)->nodeid);
//             graphno += 1;
//          }
//       }

//       void AddEdge(int refstart_i, int refend_i, int head, int tail) {
//          if (refstart_i == refstart && refend_i == refend) {
//             // no2gnode[head];
//             // fprintf(stderr, "no2gnode[0]->nodeid %d: ", no2gnode->Get(head));
//             no2gnode->Get(head)->child.Add(tail);
//             no2gnode->Get(tail)->parent.Add(head);

//             // GBitVec childpat;
//             // GBitVec parentpat;

//             edgeno+=1;
//          }
//       }

//       void PrintGraph() {
//          { //DEBUG ONLY
//             fprintf(stderr,"\tafter traverse:\n");
//             for(int i=1;i<graphno-1;i++) {
//                fprintf(stderr,"\tNode %d with parents:",i);
//                for(int p=0;p<no2gnode->Get(i)->parent.Count();p++) fprintf(stderr," %d",no2gnode->Get(i)->parent[p]);
//                fprintf(stderr," and children:");
//                for(int c=0;c<no2gnode->Get(i)->child.Count();c++) fprintf(stderr," %d",no2gnode->Get(i)->child[c]);
            
//                fprintf(stderr," %d(%d-%d)",i,no2gnode[0][i]->start,no2gnode[0][i]->end);
//                fprintf(stderr,"\n");
//             }
//          }
//       }

//       int get_refstart() {
//          return refstart;
//       }
//       int get_refend() {
//          return refend;
//       }
//       int get_s() {
//          return s;
//       }
//       int get_g_idx() {
//          return g_idx;
//       }
//       int get_graphno() {
//          return graphno;
//       }
//       int get_edgeno() {
//          return edgeno;
//       }
//       GPVec<CGraphnode>* get_no2gnode() {
//          return no2gnode;
//       }

//    // after the graph is created, I need to do 'predict transcripts for unstranded bundles here'
//    /*****************************
// 	 ** 'CPrediction': constructor
// 	 **		predict transcripts for unstranded bundles here
// 	 *****************************/

//    // And then, I need to parse the graph here
//    /*****************************
//     ** 5. parse graph
//     **    'process_refguides' & 'process_transfrags' & 'find_transcripts' & 'free_treepat'
//     *****************************/
// };

// struct UniSpliceGraphGp {
//    protected:
//       // //  Do not need
//    	// GPVec<CBundle> bundle[3]; // all bundles on all strands: 0,1,2
//       // //  Do not need
//    	// GPVec<CBundlenode> bnode[3]; // last bnodes on all strands: 0,1,2 for each bundle : this might be the key for overalps
//       // //  Do not need
//     	// int bno[2]={0,0};
//       int refstart;
// 	   int refend;
//       int gpSize[2];
//       GVec<int> graphnoGp[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
//       GVec<int> edgenoGp[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
//       GPVec<CGraphnode>* no2gnodeGp[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
//       // GVec<CGraphinfo> *bundle2graph[2]; // should I keep the neutral strand for consistency ? -> remember not to delete it
//       // //  Do not need
//       // GPVec<CTransfrag> *transfrag[2]; // for each transfrag t on a strand s, in a graph g, transfrag[s][g][t] gives it's abundance and it's pattern
//     	// // GPVec<CMTransfrag> *mgt[2]; // merged super-transfrags

//       // CTreePat **tr2no[2]; // for each graph g, on a strand s, tr2no[s][g] keeps the tree pattern structure for quick retrieval of the index t of a tansfrag
//       // GIntHash<int> *gpos[2]; // for each graph g, on a strand s, gpos[s][g] keeps the hash between edges and positions in the bitvec associated to a pattern
//       // GVec<int> lastgpos[2];
//    public:
//       // s: strand (0 = negative strand; 1 = unknown strand; 2 = positive strand // 0(-),1(.),2(+))
//       // b: all bundles on all strands: 0,1,2
//    	UniSpliceGraphGp() { 
//          refstart = 0;
// 	      refend = 0;
//          for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
//             int s=sno/2; // adjusted strand due to ignoring neutral strand
//             gpSize[s] = 0;
//             graphnoGp[s].setCapacity(8192);
//             edgenoGp[s].setCapacity(8192);
//             no2gnodeGp[s] = new GPVec<CGraphnode>[8192];
//          }
//       }

//       ~UniSpliceGraphGp() {
//          for(int i=0;i<2;i++) {
//             gpSize[i] = 0;
//             graphnoGp[i].Clear();
//             edgenoGp[i].Clear();
//             delete [] no2gnodeGp[i];
//             // no2gnodeGp[i] = new GPVec<CGraphnode>;
//             // no2gnode[i] = NULL;
//          };
//       }

//       void SetRefStartEnd(int refstart_i, int refend_i) {
//          refstart = refstart_i;
//          refend = refend_i;
//       }

//       void AddGraph(UniSpliceGraph* uni_splice_graph) {
//             // fprintf(stderr, "* uni_splice_graph.refstart: %d \n", uni_splice_graph -> get_refstart());
//             // fprintf(stderr, "* uni_splice_graph.refend: %d \n", uni_splice_graph -> get_refend());
//             // fprintf(stderr, "* uni_splice_graph.s: %d \n", uni_splice_graph -> get_s());
//             // fprintf(stderr, "* uni_splice_graph.g_idx: %d \n", uni_splice_graph -> get_g_idx());
//             // fprintf(stderr, "* uni_splice_graph.get_graphno: %d \n", uni_splice_graph -> get_graphno());
//             // fprintf(stderr, "* uni_splice_graph.get_edgeno: %d \n", uni_splice_graph -> get_edgeno());

//             // int refstart;
//             // int refend;
//             // int s;
//             // int g_idx;

//             // int graphno;  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
//             // int edgeno;  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
//             // GPVec<CGraphnode>* no2gnode; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i
//             graphnoGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()] = uni_splice_graph->get_graphno();
//             edgenoGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()] = uni_splice_graph->get_edgeno();

//             for (int i=0; i<uni_splice_graph->get_no2gnode()->Count(); i++) {
//                no2gnodeGp[uni_splice_graph->get_s()][uni_splice_graph->get_g_idx()].Add(uni_splice_graph->get_no2gnode()->Get(i));
//             }
//             gpSize[uni_splice_graph->get_s()] = uni_splice_graph->get_g_idx()+1;
//             // fprintf(stderr, "&&& uni_splice_graph->get_g_idx()+1: %d\n", uni_splice_graph->get_g_idx()+1);
//             // fprintf(stderr, "&&& gpSize[uni_splice_graph->get_s()]: %d\n", gpSize[uni_splice_graph->get_s()]);
//             // We need to reset graph_idx when (1)the new strands     
//       }

//       void Clear() {
//          // fprintf(stderr, "**** Start Clearing !!!! \n ");
//          for(int i=0;i<2;i++) {
//             gpSize[i] = 0;
//             graphnoGp[i].Clear();
//             graphnoGp[i].setCapacity(8192);
//             edgenoGp[i].Clear();
//             edgenoGp[i].setCapacity(8192);
//             delete [] no2gnodeGp[i];
//             no2gnodeGp[i] = new GPVec<CGraphnode>[8192];
//             // no2gnodeGp[i]->setCapacity(8192);
//          };
// 	   }

//       void PrintGraphGp() {
//          fprintf(stderr, "*********************************\n");
//          fprintf(stderr, "*********** PrintGraphGp ********\n");
//          fprintf(stderr, "*********************************\n");

//          { // DEBUG ONLY
//             printTime(stderr);
//             for(int s=0;s<2;s++) {
//                fprintf(stderr, "\n\tThere are %d stranded[%d] graphs\n", gpSize[s],int(2*s));
//                for(int b=0;b<gpSize[s];b++) {
//                   if(graphnoGp[s][b]) {
//                      GStr pat;
//                      fprintf(stderr,"\t\tGraph[%d][%d] with %d nodes and %d edges :",int(2*s),b,graphnoGp[s][b],edgenoGp[s][b]);
//                      for(int nd=1;nd<graphnoGp[s][b]-1;nd++)
//                         fprintf(stderr," %d(%d-%d)",nd,no2gnodeGp[s][b][nd]->start,no2gnodeGp[s][b][nd]->end);
//                   }
//                   fprintf(stderr,"\n");
//                }
//             }
//          }
//       }

//       int get_refstart () {
//          return refstart;
//       }
//       int get_refend () {
//          return refend;
//       }
//       GVec<int>* get_graphnoGp () {
//          return graphnoGp;
//       }
//       GVec<int>* get_edgenoGp () {
//          return edgenoGp;
//       }
//       GPVec<CGraphnode>** get_no2gnodeGp () {
//          return no2gnodeGp;
//       }
//       int* get_gpSize() {
//          return gpSize;
//       }
// };








class DOTReader;
class DOTWriter;

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
      
      //the caller has to FREE the created UnispgGp
      UnispgGp* next() {
         // Three target Vector to be created.
         // GVec<int> graphno[2];  // how many nodes are in a certain graph g, on strand s: graphno[s][g]
         // GVec<int> edgeno[2];  // how many edges are in a certain graph g, on strand s: edgeno[s][g]
         // GPVec<CGraphnode> *no2gnode[2]; // for each graph g, on a strand s, no2gnode[s][g][i] gives the node i



         // Example from: https://www.codegrepper.com/code-examples/cpp/remove+all+spaces+from+string+c%2B%2B
         string line;  
         if (getline(ifile_dot, line)) {
            // line.erase(remove(line.begin(), line.end(), ' '), line.end());
            // fprintf(stderr, "Before parsing line: %s \n", line.c_str());
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
            bool add_source = true;
            bool add_sink = true;
            while ((pos = line.find(delimiter)) != string::npos) {
               token = line.substr(0, pos);
               // std::cout << token << std::endl;
               // fprintf(stderr, "token: %s\n", token.c_str());
               regex title_rgx("strict\\s+digraph\\s+(\\w+)_(\\w+)_(\\w+)_(\\w+)");
               smatch match;
               fprintf(stderr, "line: %s \n", line.c_str());
               if (regex_search(token, match, title_rgx)) {
                  refstart = stoi(match[1]);
                  refend = stoi(match[2]);
                  s = stoi(match[3]);
                  g_idx = stoi(match[4]);
                  // Set back to 0 first.
                  g_idx = 0;
                  // fprintf(stderr, "ref : %d - %d, %d, %d\n", refstart, refend, s, g_idx);
                  // fprintf(stderr, "refstart: %d\n", refend);
                  // fprintf(stderr, "s: %d\n", s);
                  // fprintf(stderr, "g_idx: %d\n", g_idx);
               }

               line.erase(0, pos + delimiter.length());
               break;
            }

            // UnispgGp* uni_splice_graph = new UnispgGp();
            fprintf(stderr, ">> boundaries: %u - %u\n ", refstart, refend);
            UnispgGp* uni_splice_graph = new UnispgGp(refstart, refend);

            uni_splice_graph->ProcessSample(fname);
			   // uni_splice_graph->PrintGraphGp();

            // UniSpliceGraph* uni_splice_graph = new UniSpliceGraph(refstart, refend, s, g_idx);

            delimiter = ";";
            pos = 0;
            token = "";
            uni_splice_graph->graphno_unispg[s] = 2;
            uni_splice_graph->edgeno_unispg[s] = 0;
            while ((pos = line.find(delimiter)) != string::npos) {
               int node_id = 0;
               int start = 0;
               int end = 0;
               float cov = 0;
               int head = 0;
               int tail = 0;
               token = line.substr(0, pos);
               // std::cout << token << std::endl;
               // fprintf(stderr, "token: %s\n", token.c_str());

               regex node_rgx("(\\w+)\\[+start=(\\w+)\\s+end=(\\w+)\\s+cov=(\\w+)");
               regex edge_rgx("(\\w+)->(\\w+)");
               smatch match;

               if (add_source) {
                  // Add Source
                  GVec<bool>* is_passed_s_source = new GVec<bool>(sample_num-1, false);
                  GVec<float>* cov_s_source = new GVec<float>(sample_num-1, 0.0f);
                  GVec<float>* capacity_s_source = new GVec<float>(sample_num-1, 0.0f);
                  // delete source_gp[s];
                  CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, is_passed_s_source, cov_s_source, capacity_s_source, true, 0, 0, 0);
                  uni_splice_graph->no2gnode_unispg[s][g_idx].Add(source);
                  add_source = false;
               }

               if (regex_search(token, match, node_rgx)) {
                  node_id = stoi(match[1]);
                  start = stoi(match[2]);
                  end = stoi(match[3]);
                  cov = stof(match[4]);
                  // fprintf(stderr, "node : %d\n", node);
                  // fprintf(stderr, "start : %d\n", start);
                  // fprintf(stderr, "end : %d\n", end);
                  // fprintf(stderr, "cov : %d\n", cov);



                  // uni_splice_graph->AddNode(refstart, refend, start, end, node_id, cov);
                  // Add Node
                  GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                  GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                  GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                  // delete source_gp[s];
                  CGraphnodeUnispg* new_node = new CGraphnodeUnispg(sample_num, start, end, node_id, is_passed_s, cov_s, capacity_s, true, cov, 0, 0);
                  uni_splice_graph->no2gnode_unispg[s][g_idx].Add(new_node);
                  uni_splice_graph->graphno_unispg[s]++;
               }
               if (regex_search(token, match, edge_rgx)) {
                  // Add Sink
                  if (add_sink) {
                     GVec<bool>* is_passed_s_sink = new GVec<bool>(sample_num-1, false);
                     GVec<float>* cov_s_sink = new GVec<float>(sample_num-1, 0.0f);
                     GVec<float>* capacity_s_sink = new GVec<float>(sample_num-1, 0.0f);
                     // delete source_gp[s];
                     CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, 0, is_passed_s_sink, cov_s_sink, capacity_s_sink, true, 0, 0, 0);
                     uni_splice_graph->no2gnode_unispg[s][g_idx].Add(sink);
                     add_sink = false;
                  }
                  tail = stoi(match[1]);
                  head = stoi(match[2]);
                  // uni_splice_graph->AddEdge(refstart, refend, tail, head);

                  uni_splice_graph->no2gnode_unispg[s][g_idx].Get(tail)->child.cAdd(head);
                  uni_splice_graph->no2gnode_unispg[s][g_idx].Get(head)->parent.cAdd(tail);

                  //   node->parent.cAdd(new_parent);
                  uni_splice_graph->edgeno_unispg[s]++;
               }
               line.erase(0, pos + delimiter.length());

               // no2gnode
            }
            // uni_splice_graph -> PrintGraph();

            return uni_splice_graph;
         } else {    
            return NULL;        
         }   
      }

      bool next(UnispgGp& rec) {
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


struct DOTInputFile {
   protected:
      UnispgGp* rec;
      // UniSpliceGraphGp* uni_splice_graphs;
   public:
      DOTReader* reader;
      GStr file; //same order
      GStr tmpfile; //all the temp files created by this
      DOTInputFile():rec(NULL), reader(), file(), tmpfile() {
         // fprintf(stderr, "&& DOTInputFile Initialization  :\n");
         // uni_splice_graphs = new UniSpliceGraphGp();
         // fprintf(stderr, "Initialization uni_splice_graphs.graph_idx %d :\n", uni_splice_graphs->get_graphno());
      }
      // void Add(const char* fn);
      int count() { return 1; }
      bool start(const char* fn); //open all files, load 1 record from each
      UnispgGp* next();
      void stop();
};


#endif