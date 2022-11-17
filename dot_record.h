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
#include "helper.h"
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

      void set_samples(GVec<GStr>& samples) {
         string line;  
         if (getline(ifile_dot, line)) {
            line.erase(0,1); // removes first character
            string delimiter = ";";
            size_t pos = 0;
            string token;
            while ((pos = line.find(delimiter)) != string::npos) {
               token = line.substr(0, pos);
               line.erase(0, pos + delimiter.length());
               GStr sample_name = token.c_str();
               samples.cAdd(sample_name);
            }
         }
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

            GStr refseq = "";
            int refstart = 0;
            int refend = 0;
            int s = 0;
            int g_idx = 0;
            int g_num = 0;
            bool add_source = true;
            bool add_sink = true;
            while ((pos = line.find(delimiter)) != string::npos) {
               token = line.substr(0, pos);
               // std::cout << token << std::endl;
               // fprintf(stderr, "token: %s\n", token.c_str());
               regex title_rgx("strict\\s+digraph\\s+(\\w+)_(\\w+)_(\\w+)_(\\w+)_(\\w+)_(\\w+)");
               smatch match;
               // fprintf(stderr, "line: %s \n", line.c_str());
               if (regex_search(token, match, title_rgx)) {
                  refseq = GStr(match[1].str().c_str());
                  refstart = stoi(match[2]);
                  refend = stoi(match[3]);
                  s = stoi(match[4]);
                  g_idx = stoi(match[5]);
                  g_num = stoi(match[6]);
                  // Set back to 0 first.
                  g_idx = 0;
                  // fprintf(stderr, "refseq: %s\n", refseq.chars());
                  // fprintf(stderr, "ref : %d - %d, %d, %d\n", refstart, refend, s, g_idx);
                  // fprintf(stderr, "refstart: %d\n", refend);
                  // fprintf(stderr, "s: %d\n", s);
                  // fprintf(stderr, "g_idx: %d\n", g_idx);
               }

               line.erase(0, pos + delimiter.length());
               break;
            }

            // UnispgGp_APPLY* uni_splice_graph = new UnispgGp_APPLY();
            // fprintf(stderr, ">> boundaries: %u - %u\n ", refstart, refend);
            UnispgGp_APPLY* uni_splice_graph = new UnispgGp_APPLY(refstart, refend, g_num, refseq);
            uni_splice_graph->s_single_dot = s;

            // fprintf(stderr, ">> fname: %s\n", fname);
            // uni_splice_graph->ProcessSample(fname);
			   // uni_splice_graph->PrintGraphGp();

            // UniSpliceGraph* uni_splice_graph = new UniSpliceGraph(refstart, refend, s, g_idx);

            delimiter = ";";
            pos = 0;
            token = "";
            int node_num = 0;
            int edge_num = 0;
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
                  node_num += 1;
                  // Add Source
                  // GVec<bool>* is_passed_s_source = new GVec<bool>(sample_num-1, false);
                  // GVec<float>* cov_s_source = new GVec<float>(sample_num-1, 0.0f);
                  // GVec<float>* capacity_s_source = new GVec<float>(sample_num-1, 0.0f);

                  // delete source_gp[s];
                  CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, 0, NULL, NULL, NULL, true, 0, 0, 0);
                  uni_splice_graph->no2gnode_unispg[s][g_idx].Add(source);
                  add_source = false;
               }

               if (regex_search(token, match, node_rgx)) {
                  node_id = stoi(match[1]);
                  // if (max_node_id < node_id) max_node_id = node_id;
                  node_num += 1;
                  start = stoi(match[2]);
                  end = stoi(match[3]);
                  cov = stof(match[4]);
                  // fprintf(stderr, "node : %d\n", node);
                  // fprintf(stderr, "start : %d\n", start);
                  // fprintf(stderr, "end : %d\n", end);
                  // fprintf(stderr, "cov : %d\n", cov);



                  // uni_splice_graph->AddNode(refstart, refend, start, end, node_id, cov);
                  // Add Node
                  // GVec<bool>* is_passed_s = new GVec<bool>(sample_num-1, false);
                  // GVec<float>* cov_s = new GVec<float>(sample_num-1, 0.0f);
                  // GVec<float>* capacity_s = new GVec<float>(sample_num-1, 0.0f);
                  // delete source_gp[s];
                  CGraphnodeUnispg* new_node = new CGraphnodeUnispg(sample_num, start, end, node_id, NULL, NULL, NULL, true, cov, 0, 0);
                  uni_splice_graph->no2gnode_unispg[s][g_idx].Add(new_node);
               }
               if (regex_search(token, match, edge_rgx)) {
                  // Add Sink
                  if (add_sink) {
                     node_num += 1;
                     // GVec<bool>* is_passed_s_sink = new GVec<bool>(sample_num-1, false);
                     // GVec<float>* cov_s_sink = new GVec<float>(sample_num-1, 0.0f);
                     // GVec<float>* capacity_s_sink = new GVec<float>(sample_num-1, 0.0f);
                     // delete source_gp[s];
                     CGraphnodeUnispg* sink = new CGraphnodeUnispg(sample_num, 0, 0, node_num-1, NULL, NULL, NULL, true, 0, 0, 0);
                     uni_splice_graph->no2gnode_unispg[s][g_idx].Add(sink);
                     add_sink = false;
                  }
                  edge_num += 1;

                  tail = stoi(match[1]);
                  head = stoi(match[2]);

                  // Consturcting edge key.
                  // fprintf(stderr, ">>> node_num during construction: %d;  edge_num: %d\n", node_num, edge_num);
				      int key=edge(tail, head, node_num);
                  int *pos=uni_splice_graph->gpos[s][g_idx][key];

						// fprintf(stderr, ">>> @@ Adding edge (%d - %d);  g_node_num: %d\n", tail, head, node_num);
						// fprintf(stderr, ">>> *pos: %d\n", *pos);

                  if(pos==NULL) {
                     uni_splice_graph->gpos[s][g_idx].Add(key, node_num-1+edge_num);
						   // fprintf(stderr, ">>> (1) key: %d; *pos: %d\n", key, node_num-1+edge_num);
                  } else {
						   // fprintf(stderr, ">>> (2) key: %d; *pos: %d\n", key, *pos);
                  }
                  // uni_splice_graph->AddEdge(refstart, refend, tail, head);

                  uni_splice_graph->no2gnode_unispg[s][g_idx].Get(tail)->child.cAdd(head);
                  uni_splice_graph->no2gnode_unispg[s][g_idx].Get(head)->parent.cAdd(tail);

                  //   node->parent.cAdd(new_parent);
               }
               line.erase(0, pos + delimiter.length());

               // no2gnode
            }
            // uni_splice_graph -> PrintGraph();
            uni_splice_graph->set_N_E_num(s, node_num, edge_num);
            uni_splice_graph->graph_num[s] += 1;
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
		GVec<GStr> samples;
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
      void set_samples();
      UnispgGp_APPLY* next();
      void stop();
};

#endif