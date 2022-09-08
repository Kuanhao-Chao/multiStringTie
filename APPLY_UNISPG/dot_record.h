#ifndef __DOT_RECORD_H__
#define __DOT_RECORD_H__
// #include "rlink.h"
// #include "GArgs.h"
// #include "GStr.h"

// #include "gff.h"
// #include "GSam.h"
// #include "GBitVec.h"
// #include "time.h"
// #include "tablemaker.h"
// #include "GHashMap.hh"
// #include "unispg.h"
#include "global_params.h"
#include "APPLY_UNISPG/unispg_A.h"

#include <fstream>
#include <regex>
#include <string>
#include <typeinfo>

using namespace std;

// class UnispgGp_APPLY;
// class DOTReader;
// class DOTWriter;

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
      
      //the caller has to FREE the created UnispgGp_APPLY
      UnispgGp_APPLY* next() {
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

            // UnispgGp_APPLY* uni_splice_graph = new UnispgGp_APPLY();
            fprintf(stderr, ">> boundaries: %u - %u\n ", refstart, refend);
            UnispgGp_APPLY* uni_splice_graph = new UnispgGp_APPLY(refstart, refend);

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

                  GVec<float>* capacity_s_unispg_source = new GVec<float>(sample_num-1, 0.0f);

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

      bool next(UnispgGp_APPLY& rec) {
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
      UnispgGp_APPLY* rec;
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
      UnispgGp_APPLY* next();
      void stop();
};

#endif