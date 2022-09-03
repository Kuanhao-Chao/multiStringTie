#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "GHashMap.hh"

#include "tmerge.h"
#include "multist.h"
#include "unispg.h"
#include "helper.h"
#include "processOptions.h"
#include "processBundle.h"
#include "definitions.h"

// Include fstream to quickly write out the ratio!!!
#include <fstream>
#include <vector> 
#include "time.h"
#include <iostream>

/*******************************************
 ** Argument parsing parameters.
 *******************************************/
bool debugMode=false; // "debug" or "D" tag.
bool verbose=false; // "v" tag.
bool keepTempFiles; // "keeptmp" tag.
bool fr_strand=false; // "fr" tag.
bool rf_strand=false; // "rf" tag.
bool guided=false; // "G" tag.
bool viral=false; // "viral" tag.
bool eonly=false; // "e" tag. for mergeMode includes estimated coverage sum in the merged transcripts
bool longreads=false; // "L" tag.
GStrSet<> excludeGseqs; // "x" tag. hash of chromosomes/contigs to exclude (e.g. chrM)
uint bundledist=50;  // "g" tag. reads at what distance should be considered part of separate bundles
uint runoffdist=200; // threshold for 'bundledist'
float readthr=1; // "c" tag. read coverage per bundle bp to accept it; // paper uses 3
float tpm_thr=1; // "T" tag.
float fpkm_thr=1; // "F" tag.
bool isunitig=true; // "U" tag.
uint junctionsupport=10; // "a" tag. anchor length for junction to be considered well supported <- consider shorter??
uint sserror=25; // "E" tag. window arround splice sites that we use to generate consensus in case of long read data
int junctionthr=1; // "j" tag. number of reads needed to support a particular junction
float mcov=1; // "M" tag. fraction of bundle allowed to be covered by multi-hit reads paper uses 1
int mintranscriptlen=200; // "m" tag. minimum length for a transcript to be printed
float singlethr=4.75; // "s" tag. coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
bool mixedMode=false; // "mix" tag. both short and long read data alignments are provided
float isofrac=0.01; // "f" tag. 
bool trim=true; // "t" tag. 
bool includesource=true; // "z" tag.
bool nomulti=false; // "u" tag.
bool multiMode=false; // "multi" tag.
bool mergeMode = false; // "merge" tag. 
GFastaDb* gfasta=NULL; // "rseq" or "S" tag.
int num_cpus=1; // "p" tag.
bool rawreads=false; // "R" tag.
GStr label("STRG"); // "l" tag.
bool retained_intron=false; // "i" tag. set by parameter -i for merge option
bool geneabundance=false; // "A" tag.
GStr guidegff; // "G" tag
GStr ptff; // "ptf" tag. (point features)
FILE* f_out=NULL; // get from "outfname" param.
FILE* c_out=NULL; // "C" tag.


/*******************************************
 ** File-related parameters.
 *******************************************/
// Dot file => for reading in DOT format.
GStr unispgdotfname_root; 
GStr unispgdotfname_pos; 
GStr unispgdotfname_neg; 
// Annotation file => outputing in GFF format.
GStr outfname;
GStr out_dir;
GStr tmp_path;
GStr cram_ref; //reference genome FASTA for CRAM input
GStr tmpfname;
GStr genefname;
GStr traindir; // training directory for CDS option
// Ratio file => for coverage comparison & visualization.
ofstream ratio_file_pos;
ofstream ratio_file_neg;

/*******************************************
 ** Program-defined global parameter.
 *******************************************/
multiStringTieMode mode;
bool NoMoreBundles=false;
GffNames* gseqNames=NULL; //used as a dictionary for reference sequence names and ids
int refseqCount=0; // number of reference sequences found in the guides file
bool includecov=false;
// Global counter variable.
int sample_num = 1;
double Num_Fragments=0; //global fragment counter (aligned pairs)
double Frag_Len=0;
double Cov_Sum=0;
int GeneNo=0; //-- global "gene" counter
// For multistringtie applyUNISPG only. Counting how many reads are skipped (not processed).
int skip_counter = 0;

/*******************************************
 ** Reader parameters.
 *******************************************/
TInputFiles bamreader;
DOTInputFile dotreader;
DOTInputFile dotreader_pos;
DOTInputFile dotreader_neg;





#ifndef NOTHREADS
//single producer, multiple consumers
//main thread/program is always loading the producer
GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
int bundleWork=1; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue

//GFastMutex waitMutex;
GMutex waitMutex; // controls threadsWaiting (idle threads counter)

int threadsWaiting; // idle worker threads
GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle

GConditionVar haveClear; //will notify when bundle buf space available

GMutex queueMutex; //controls bundleQueue and bundles access

GFastMutex printMutex; //for writing the output to file

GFastMutex logMutex; //only when verbose - to avoid mangling the log output

GFastMutex bamReadingMutex;

GFastMutex countMutex;

#endif

int main(int argc, char*argv[]) {
	std::cout << "This is the multiStringTie program.\n" << std::endl;

	fprintf(stderr, "%d  %s  %s  %s  %s\n", argc, argv[0], argv[1], argv[2], argv[3]);
	if (strcmp(argv[1], "CREATE_UNISPG") || strcmp(argv[1], "APPLY_UNISPG")) {
		// fprintf(stderr, "This is mode: %s  %d\n", argv[1], strcmp(argv[1], "APPLY_UNISPG"));
		if (strcmp(argv[1], "CREATE_UNISPG") == 0) {
			mode = CREATE_UNISPG;
		}
		if (strcmp(argv[1], "APPLY_UNISPG") == 0) {
			mode = APPLY_UNISPG;
		}
		for (int i=1; i<argc-1; i++) {
			argv[i] = argv[i+1];
		}
		argv[argc-1] = NULL;
		argc = argc-1;
		fprintf(stderr, "%d  %s  %s  %s  %s\n", argc, argv[0], argv[1], argv[2], argv[3]);
	}
   
   
	if (mode == CREATE_UNISPG) {
		GArgs args(argc, argv,
		"debug;help;version;viral;conservative;mix;unispg=;ref=;cram-ref=cds=;keeptmp;rseq=;ptf=;bam;fr;rf;merge;multi;"
		"exclude=zihvteuLRx:n:j:s:D:G:C:S:l:m:o:a:j:c:f:p:g:P:M:Bb:A:E:F:T:");
		args.printError(USAGE, true);
		processCreateOptions(args);


	} else if (mode == APPLY_UNISPG) {
		GArgs args(argc, argv,
		"debug;help;version;dot=;d=;bam;fr;rf;"
		"exclude=hvdx:n:j:s:d:D:G:C:S:l:m:o:a:j:c:f:p:g:P:M:Bb:A:E:F:T:");

		processApplyOptions(args);

		GVec<int> alncounts(30); //keep track of the number of read alignments per chromosome [gseq_id]

		int bamcount=bamreader.start(); //setup and open input files
		fprintf(stderr, "&& bamcount number: %d\n", bamcount);

	#ifndef GFF_DEBUG
		if (bamcount<1) {
			GError("%sError: no input BAM provided!\n", USAGE);
		} else if (bamcount > 1) {
			GError("%sError: more than 1 BAM file is provided.\n", USAGE);
		}
	#endif

	#ifdef DEBUGPRINT
		verbose=true;
	#endif
		const char* ERR_BAM_SORT="\nError: the input alignment file is not sorted!\n";



		/*******************************************
		 *******************************************
		 ** input processing
		 *******************************************
		 *******************************************/
		GHash<int> hashread;      //read_name:pos:hit_index => readlist index
		GList<GffObj>* guides=NULL; //list of transcripts on a specific reference
		GList<GPtFeature>* refptfs=NULL; //list of point-features on a specific reference
		// int currentstart=0, currentend=0;
		int ng_start=0;
		int ng_end=-1;
		int ptf_idx=0; //point-feature current index in the current (*refptfs)[]
		int ng=0;
		GStr lastref;
		bool no_ref_used=true;
		int lastref_id=-1; //last seen gseq_id
		// int ncluster=0; used it for debug purposes only

	#ifndef NOTHREADS

		//model: one producer, multiple consumers
	#define DEF_TSTACK_SIZE 8388608
		size_t defStackSize=DEF_TSTACK_SIZE;
	#ifdef _GTHREADS_POSIX_
		int tstackSize=GThread::defaultStackSize();
		if (tstackSize<DEF_TSTACK_SIZE) defStackSize=DEF_TSTACK_SIZE;
		if (verbose) {
		if (defStackSize>0){
			int ssize=defStackSize;
			GMessage("Default stack size for threads: %d (increased to %d)\n", tstackSize, ssize);
		}
		else GMessage("Default stack size for threads: %d\n", tstackSize);
		}
	#endif
		GThread* threads=new GThread[num_cpus]; //bundle processing threads

		GPVec<BundleData> bundleQueue(false); //queue of loaded bundles
		//the consumers take (pop) bundles out of this queue for processing
		//the producer populates this queue with bundles built from reading the BAM input

		BundleData* bundles=new BundleData[num_cpus+1];
		//bundles[0..num_cpus-1] are processed by threads, loading bundles[num_cpus] first

		dataClear.setCapacity(num_cpus+1);
		for (int b=0;b<num_cpus;b++) {
			threads[b].kickStart(workerThread, (void*) &bundleQueue, defStackSize);
			bundles[b+1].idx=b+1;
			dataClear.Push(b);
		}
		BundleData* bundle = &(bundles[num_cpus]);
	#else
		BundleData bundles[1];
		BundleData* bundle = &(bundles[0]);

		// if (multiMode) {
		// } else {
		// 	UniSpliceGraphGp uni_splice_graphGps[1];
		// 	UniSpliceGraphGp* uni_splice_graphGp = &(uni_splice_graphGps[0]);
		// }
		// UniSpliceGraphGp uni_splice_graphGps[1];
		// UniSpliceGraphGp* uni_splice_graphGp = &(uni_splice_graphGps[0]);
	#endif

		/*******************************************
		 ** This is for reading bam file.
		 *******************************************/
		GSamRecord* brec=NULL;	
		GSamRecord* prev_brec=NULL;	
		bool more_alns=true;
		TAlnInfo* tinfo=NULL; // for --merge

		/*******************************************
		 ** This is for reading dot file.
		 *******************************************/
		bool more_graph=true;
		bool more_graph_pos=true;
		bool more_graph_neg=true;

		bool next_graph_pos=true;
		bool next_graph_neg=true;

		int prev_pos=0;
		// bool skipGseq=false;
		UnispgGp* drec=NULL;
		UnispgGp* prev_drec=NULL;
		UnispgGp* drec_pos=NULL;
		UnispgGp* drec_neg=NULL;
		// UnispgGp* prev_drec_pos=NULL;
		// UnispgGp* prev_drec_neg=NULL;
		int prev_drec_pos_start=NULL;
		int prev_drec_pos_end=NULL;
		int prev_drec_neg_start=NULL;
		int prev_drec_neg_end=NULL;

		bool next_pos_neg = NULL;

		// Initialize the dot graphs vector.
		GPVec<UnispgGp>* graphs_vec[2];
		for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
			int s=sno/2; // adjusted strand due to ignoring neutral strand
			graphs_vec[s] = new GPVec<UnispgGp>[1024];
		}


		// while ((brec=bamreader.next())!=NULL) {
		// }




		fprintf(stderr, "unispgdotfname_root: %s\n", unispgdotfname_root.chars());

		unispgdotfname_pos = unispgdotfname_root + "_node_pos_0_unispg.dot";
		unispgdotfname_neg = unispgdotfname_root + "_node_neg_0_unispg.dot";
		fprintf(stderr, "unispgdotfname_pos: %s\n", unispgdotfname_pos.chars());
		fprintf(stderr, "unispgdotfname_neg: %s\n", unispgdotfname_neg.chars());

		// unispgdotfname_root
		if (fileExists(unispgdotfname_pos.chars())==0 || fileExists(unispgdotfname_neg.chars())==0) {
			GError("%sError: the dot file does not exist!\n", USAGE);
		}

		bool dot_is_open = dotreader.start(unispgdotfname_pos); //setup and open DOT input file
		fprintf(stderr, "unispgdotfname_pos: %d\n", dot_is_open);


		bool dot_is_open_pos = dotreader_pos.start(unispgdotfname_pos); //setup and open DOT input file
		fprintf(stderr, "unispgdotfname_pos: %d\n", dot_is_open_pos);

		bool dot_is_open_neg = dotreader_neg.start(unispgdotfname_neg); //setup and open DOT input file
		fprintf(stderr, "unispgdotfname_neg: %d\n", dot_is_open_neg);

		fprintf(stderr, "bamreader.files.Get(file_idx): %s\n", bamreader.files.Get(0).chars());


		bundle->Clear();
		brec=NULL;	
		more_alns=true;
		more_graph=true;
		more_graph_pos=true;
		more_graph_neg=true;
		next_graph_pos=true;
		next_graph_neg=true;

		tinfo=NULL; // for --merge
		prev_pos=0;
		// skipGseq=false;
		uint r_start = 0;
		uint r_end = 0;
		char xstrand='.';
		const char* refseqName=NULL;
		int nh=1;
		int hi=0;
		int gseq_id=lastref_id;  //current chr id
		bool new_bundle=false;


		fprintf(stderr, "** file_idx: %d\n", 0);
		bamreader.start_fidx(0);



		int total_read = 0;
		int processed_read = 0;
		int pos_strand = 0;
		int neg_strand = 0;
		int unstrand = 0;

		bool process_graphs_vec = false;


  		ratio_file_pos.open("/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/ratio_pos.txt");
  		ratio_file_neg.open("/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/ratio_neg.txt");



		/*****************************
		 * Processing Graphs & alignments: main algorithm
		 *****************************/
		/*****************************
		 * Iterating graphs
		 *****************************/
		while (more_graph_pos || more_graph_neg) {
			// After this, the current drec_pos is the newest positive graph.
			if (next_graph_pos) {
				if ((drec_pos=dotreader_pos.next())!=NULL) {
					fprintf(stderr, "***** drec_pos: %d - %d\n", drec_pos->get_refstart(), drec_pos->get_refend());
				} else {
					more_graph_pos = false;
					next_graph_pos = false;
				}
			}
			// After this, the current drec_pos is the newest negative graph.
			if (next_graph_neg) {
				if ((drec_neg=dotreader_neg.next())!=NULL) {
					fprintf(stderr, "***** drec_neg: %d - %d\n", drec_neg->get_refstart(), drec_neg->get_refend());
				} else {
					more_graph_neg = false;
					next_graph_neg = false;
				}
			}

			if ( (prev_drec_neg_start == NULL && prev_drec_neg_end == NULL && prev_drec_pos_start == NULL && prev_drec_pos_end == NULL) || (graphs_vec[0]->Count()==0 && graphs_vec[1]->Count()==0) ) {
				// This is for beginning condition.
				fprintf(stderr, ">> This is for beginning condition.\n");
				fprintf(stderr, ">> Both `prev_drec_neg` and `prev_drec_pos` are NULL.\n");
				if (more_graph_pos && more_graph_neg) {
					fprintf(stderr, ">> drec_pos->get_refstart(): %u \n>> drec_neg->get_refstart(): %u \n", drec_pos->get_refstart(), drec_neg->get_refstart());
					process_graphs_vec = false;
					if (drec_pos->get_refstart() <= drec_neg->get_refstart()) {
						fprintf(stderr, "Adding `drec_pos`!!!!\n");

						// UnispgGp* copy_drec_pos = new UnispgGp(drec_pos);
						graphs_vec[1]->Add(drec_pos);
						prev_drec_pos_start = drec_pos->get_refstart();
						prev_drec_pos_end = drec_pos->get_refend();

						next_graph_pos = true;
						next_graph_neg = false;
					} else if (drec_pos->get_refstart() > drec_neg->get_refstart()) {
						fprintf(stderr, "Adding `drec_neg`!!!!\n");
						
						// UnispgGp* copy_drec_neg = new UnispgGp(drec_neg);
						graphs_vec[0]->Add(drec_neg);
						prev_drec_neg_start = drec_neg->get_refstart();
						prev_drec_neg_end = drec_neg->get_refend();

						next_graph_pos = false;
						next_graph_neg = true;
					}
				} else {
					process_graphs_vec = false;
					break;
				}
			} else {
				// graphs_vec is not empty in the beginning => start chaning now!
				if (prev_drec_neg_start != NULL && prev_drec_neg_end != NULL) {
					fprintf(stderr, ">> pre_neg, cur_pos\n");
					fprintf(stderr, ">> prev_drec_neg_start: %d; prev_drec_neg_end: %d\n", prev_drec_neg_start, prev_drec_neg_end);
					if (more_graph_pos) {
						bool overlap = segs_overlap(prev_drec_neg_start, prev_drec_neg_end, drec_pos->get_refstart(), drec_pos->get_refend());
						fprintf(stderr, ">> overlap: %d", overlap);

						if (overlap) {
							process_graphs_vec = false;
							next_graph_pos = true;
							fprintf(stderr, ">> Adding `drec_pos`: %d - %d \n", drec_pos->get_refstart(), drec_pos->get_refend());
							fprintf(stderr, "\t\t>> `prev_drec_neg`: %d - %d \n", prev_drec_neg_start, prev_drec_neg_end);
							// UnispgGp* copy_drec_pos = new UnispgGp(drec_pos);
							graphs_vec[1]->Add(drec_pos);
							prev_drec_pos_start = drec_pos->get_refstart();
							prev_drec_pos_end = drec_pos->get_refend();
							next_graph_pos = true;
							next_graph_neg = false;
							continue;
						} else {
							process_graphs_vec = true;
							next_graph_pos = false;
							next_graph_neg = false;
						}
					}
					fprintf(stderr, ">> process_graphs_vec: %d", process_graphs_vec);
				}

				if (prev_drec_pos_start != NULL && prev_drec_pos_end != NULL) {
					fprintf(stderr, ">> pre_pos, cur_neg\n");
					if (more_graph_neg) {
						bool overlap = segs_overlap(prev_drec_pos_start, prev_drec_pos_end, drec_neg->get_refstart(), drec_neg->get_refend());
						if (overlap) {
							process_graphs_vec = false;
							next_graph_neg = true;
							fprintf(stderr, ">> Adding `drec_neg`: %d - %d \n", drec_neg->get_refstart(), drec_neg->get_refend());
							fprintf(stderr, "\t\t>> `prev_drec_pos`: %d - %d \n", prev_drec_pos_start, prev_drec_pos_end);
							// UnispgGp* copy_drec_neg = new UnispgGp(drec_neg);
							graphs_vec[0]->Add(drec_neg);
							prev_drec_neg_start = drec_neg->get_refstart();
							prev_drec_neg_end = drec_neg->get_refend();
							next_graph_pos = false;
							next_graph_neg = true;
							continue;
						} else {
							process_graphs_vec = true;
							next_graph_pos = false;
							next_graph_neg = false;
						}
					}
				}
			}



			/*****************************
			 * The graphs_vec overlapping checking is done. There are no overlaps between pos / neg.  Now, we need to process it.
			 *****************************/
			if (process_graphs_vec) {
				int graphs_vec_start = 0;
				int graphs_vec_end = 0;

				if (graphs_vec[0]->First() == nullptr) {
					graphs_vec_start = graphs_vec[1]->First()->refstart;
				} else if (graphs_vec[1]->First() == nullptr) {
					graphs_vec_start = graphs_vec[0]->First()->refstart;
				} else {
					graphs_vec_start = (graphs_vec[0]->First()->refstart < graphs_vec[1]->First()->refstart) ? graphs_vec[0]->First()->refstart : graphs_vec[1]->First()->refstart;
				}

				if (graphs_vec[0]->Last() == nullptr) {
					graphs_vec_end = graphs_vec[1]->Last()->refend;
				} else if (graphs_vec[1]->Last() == nullptr) {
					graphs_vec_end = graphs_vec[0]->Last()->refend;
				} else {
					graphs_vec_end = (graphs_vec[0]->Last()->refend > graphs_vec[1]->Last()->refend) ? graphs_vec[0]->Last()->refend : graphs_vec[1]->Last()->refend;
				}

				fprintf(stderr, "\tgraphs_vec_boundary: %d - %d \n", graphs_vec_start, graphs_vec_end);


				/*****************************
				 * Reading BAM file & check reads overlapping with the `graphs_vec`
				 *****************************/
				bool read_in_unispg = true;
				while (more_alns && read_in_unispg) {
					total_read += 1;
					bool chr_changed=false;
					int pos=0;
					const char* refseqName=NULL;
					char xstrand=0;
					int nh=1;
					int hi=0;
					int gseq_id=lastref_id;  //current chr id
					bool new_bundle=false;
					bool process_read=true;
					//delete brec;
					if ((brec=bamreader.next())!=NULL) {
						if (brec->isUnmapped()) continue;
						if (brec->start<1 || brec->mapped_len<10) {
							if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
									brec->name(), brec->start, brec->mapped_len);
							continue;
						}
#ifdef DBG_ALN_DATA
						dbg_waln(brec);
#endif

						r_start = brec->start; //start<end always!
						r_end = brec->end;

						refseqName=brec->refName();
						xstrand=brec->spliceStrand(); // tagged strand gets priority

						// fprintf(stderr, "** refseqName：%s\n", refseqName);
						processed_read += 1;
						/*****************************
						 * set strand if stranded library
						 *****************************/
						fprintf(stderr, "** Before setting strand: %c\n", xstrand);
						if(xstrand=='.' && (fr_strand || rf_strand)) {
							if(brec->isPaired()) { // read is paired
								if(brec->pairOrder()==1) { // first read in pair
									if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
									else xstrand='-';
								}
								else {
									if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
									else xstrand='+';
								}
							}
							else {
								if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
								else xstrand='-';
							}
						}
						fprintf(stderr, "** After setting strand: %c\n", xstrand);


						/*****************************
						 * Found the overlapping between a read & the unispg
						 *****************************/
						// if (xstrand == '+' || xstrand=='.') {
						if (xstrand == '+') {
							pos_strand += 1;
						} else if (xstrand=='-') {
							neg_strand += 1;
						} else if (xstrand=='.') {
							unstrand += 1;
							fprintf(stderr, "Add new unstrand read!!\n");
						}
						/*****************************
						 ** Step 1-2: Check whether reads are in the range of the graph.
						*****************************/
						// fprintf(stderr, "Process read (pre: %d - %d ;  now: %d - %d)!!!\n", pre_refstart, pre_refend, brec->start, brec->end);
						fprintf(stderr, "\t>>>>> brec: %d - %d\n", brec->start, brec->end);
						if (brec->end <  graphs_vec_start) {
							// ----------   |(s).................(e)|
							// The read is outside the current bundle => skipped!
							fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
							process_read = false;
							// continue;
							// read_in_unispg = false;
							// new_bundle = true;
						}
						if (brec->start < graphs_vec_start && brec->end == graphs_vec_start) {
							// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
							fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
							// continue;
						} else if (brec->start < graphs_vec_start && brec->end > graphs_vec_start) {
							// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
							fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
						} else if (brec->start < graphs_vec_start && brec->end == graphs_vec_start) {
							// -----|(s)---------(e)|
							fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
						} else if (brec->start == graphs_vec_start && brec->start < graphs_vec_end && brec->end < graphs_vec_end) {
							// |(s)----------.................(e)|   or   |(s)....----------........(e)|
							fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
						} else if (brec->start == graphs_vec_start && brec->start < graphs_vec_end && brec->end == graphs_vec_end) {
							// |(s)----------.................(e)|   or   |(s)....----------........(e)|
							fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
						} else if (brec->start == graphs_vec_start && brec->start < graphs_vec_end && brec->end > graphs_vec_end) {
							// |(s)----------.................(e)|   or   |(s)....----------........(e)|
							fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
						} else if (brec->start > graphs_vec_start && brec->start < graphs_vec_end && brec->end < graphs_vec_end) {
							// |(s)----------.................(e)|   or   |(s)....----------........(e)|
							fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
						} else if (brec->start > graphs_vec_start && brec->start < graphs_vec_end && brec->end == graphs_vec_end) {
							// |(s)----------.................(e)|   or   |(s)....----------........(e)|
							fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
						} else if (brec->start > graphs_vec_start && brec->start < graphs_vec_end && brec->end > graphs_vec_end) {
							// |(s)----------.................(e)|   or   |(s)....----------........(e)|
							fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
						} else if (brec->start > graphs_vec_start && brec->start < graphs_vec_end && brec->end > graphs_vec_end) {
							// |(s)...............------(e)|-----    or   |(s).................(e)|----------   
							// The overlapping with the current processing bundle.
							fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
						} else if (brec->start > graphs_vec_start && brec->start == graphs_vec_end && brec->end > graphs_vec_end) {
							// |(s)...............------(e)|-----    or   |(s).................(e)|----------   
							// The overlapping with the current processing bundle.
							fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
							// continue;
							// int overlap_current = 0;
							// overlap_current = drec->get_refend() - brec->start + 1;

							// int overlap_next = 0;
							// overlap_next = brec->end - drec->get_refend() + 1;

							// if (overlap_current > overlap_next) {
							// 	// The read belongs to the current processing bundle.
							// } else {
							// 	// The read belongs to the next bundle or not belongs to any bundles.
							// 	read_in_unispg = false;
							// }
						} else if (brec->start > graphs_vec_end) {
							fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
							read_in_unispg = false;
							new_bundle = true;
							process_read = false;
						}
					} else { //no more alignments
						more_alns=false;
						new_bundle=true; //fake a new start (end of last bundle)
					}

					/*****************************
					 * Process the new bundle
					 *****************************/
					if (new_bundle) {
						fprintf(stderr, "This is a new bundle\n");
						fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
						hashread.Clear();
						if (bundle->readlist.Count()>0) { // process reads in previous bundle
							// (readthr, junctionthr, mintranscriptlen are globals)
							bundle->getReady(graphs_vec_start, graphs_vec_end);
				#ifndef NOTHREADS
							//push this in the bundle queue where it'll be picked up by the threads
							DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
							int qCount=0;
							queueMutex.lock();
							bundleQueue.Push(bundle);
							bundleWork |= 0x02; //set bit 1
							qCount=bundleQueue.Count();
							queueMutex.unlock();
							DBGPRINT2("##> bundleQueue.Count()=%d)\n", qCount);
							//wait for a thread to pop this bundle from the queue
							waitMutex.lock();
							DBGPRINT("##> waiting for a thread to become available..\n");
							while (threadsWaiting==0) {
								haveThreads.wait(waitMutex);
							}
							waitMutex.unlock();
							haveBundles.notify_one();
							DBGPRINT("##> waitMutex unlocked, haveBundles notified, current thread yielding\n");
							current_thread::yield();
							queueMutex.lock();
							DBGPRINT("##> queueMutex locked until bundleQueue.Count()==qCount\n");
							while (bundleQueue.Count()==qCount) {
								queueMutex.unlock();
								DBGPRINT2("##> queueMutex unlocked as bundleQueue.Count()==%d\n", qCount);
								haveBundles.notify_one();
								current_thread::yield();
								queueMutex.lock();
								DBGPRINT("##> queueMutex locked again within while loop\n");
							}
							queueMutex.unlock();
				#else //no threads
							//Num_Fragments+=bundle->num_fragments;
							//Frag_Len+=bundle->frag_len;
							processBundleUnispg(bundle, graphs_vec);
				#endif
							// ncluster++; used it for debug purposes only
						} //have alignments to process
						else { //no read alignments in this bundle?
				#ifndef NOTHREADS
							dataMutex.lock();
							DBGPRINT2("##> dataMutex locked for bundle #%d clearing..\n", bundle->idx);
				#endif
							bundle->Clear();
				#ifndef NOTHREADS
							dataClear.Push(bundle->idx);
							DBGPRINT2("##> dataMutex unlocking as dataClear got pushed idx #%d\n", bundle->idx);
							dataMutex.unlock();
				#endif
						} //nothing to do with this bundle
						// if (chr_changed) {
						// 	lastref=refseqName;
						// 	lastref_id=gseq_id;
						// 	// currentend=0;
						// }
						if (!more_alns) {
							if (verbose) {
				#ifndef NOTHREADS
								GLockGuard<GFastMutex> lock(logMutex);
				#endif
								if (Num_Fragments) {
									// printTime(stderr);
									GMessage(" %g aligned fragments found.\n", Num_Fragments);
								}
								//GMessage(" Done reading alignments.\n");
							}
							noMoreBundles();
							break;
						}
			#ifndef NOTHREADS
						int new_bidx=waitForData(bundles);
						if (new_bidx<0) {
							//should never happen!
							GError("Error: waitForData() returned invalid bundle index(%d)!\n",new_bidx);
							break;
						}
						bundle=&(bundles[new_bidx]);
			#endif
						// currentstart=brec->start;
						// currentend=brec->end;
						bundle->refseq=lastref;
						bundle->start=brec->start;
						bundle->end= brec->end;
					}
				
					/*****************************
					 * Actually start processing reads.
					 *****************************/
					if (brec != NULL) {
						nh=brec->tag_int("NH");
						if (nh==0) nh=1;
						hi=brec->tag_int("HI");
						// Old conditions to stop adding reads into a bundle
						// if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
						// 	new_bundle=true;
						// }
						fprintf(stderr, "brec: %d\n", brec->mapped_len);
						// if (brec->start > drec->get_refend()) {
						fprintf(stderr, "brec: %d - %d\n", brec->start, brec->end);
						// 	new_bundle=true;
						// }
						GReadAlnData alndata(brec, 0, nh, hi, tinfo);
						fprintf(stderr, ">>> before evalReadAln xstrand: %c\n", xstrand);
						bool ovlpguide=bundle->evalReadAln(alndata, xstrand);
						fprintf(stderr, ">>> after evalReadAln xstrand: %c\n", xstrand);

						if (xstrand=='+') alndata.strand=1;
						else if (xstrand=='-') alndata.strand=-1;
						//GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
						//countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
						//fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);
						// processRead(brec->start, brec->end, *bundle, hashread, alndata);
						processRead(graphs_vec_start, graphs_vec_end, *bundle, hashread, alndata);
						
						if (new_bundle) {
							break;	
						}
					}
				}

				/*****************************
				 * Clean up the 'graphs_vec'.
				 *****************************/
				for(int s=0;s<2;s+=1) { // skip neutral bundles -> those shouldn't have junctions
					// Before cleaning up, check the total number!!!
					fprintf(stderr, "\t@@ cleanup graphs_vec[%d]: %d\n", s, graphs_vec[s]->Count());
					for (int i=0; i<graphs_vec[s]->Count(); i++) {
						fprintf(stderr, "\t\t bound: %u - %u \n", graphs_vec[s]->Get(i)->get_refstart(), graphs_vec[s]->Get(i)->get_refend());
					}

					// fprintf(stderr, ">> graphs_vec[%d] Capacity: %d\n", s, graphs_vec[s]->Capacity());

					// delete [] graphs_vec[s];
					graphs_vec[s]->Clear();
					graphs_vec[s] = new GPVec<UnispgGp>[1024];
				}
			}
		}
		
		fprintf(stderr, "total read num: %d\n", total_read);
		fprintf(stderr, "processed read num: %d\n", processed_read);
		fprintf(stderr, "positve read num: %d\n", pos_strand);
		fprintf(stderr, "negative read num: %d\n", neg_strand);
		fprintf(stderr, "unstranded read num: %d\n", unstrand);
		fprintf(stderr, "skip_counter: %d\n", skip_counter);
		delete brec;
		bamreader.stop(); //close all BAM files
		delete drec;
		dotreader.stop(); //close all DOT files
		// dotreader_pos.stop(); //close all DOT files
		// dotreader_neg.stop(); //close all DOT files

  		ratio_file_pos.close();
  		ratio_file_neg.close();
	}
    return 0;
}