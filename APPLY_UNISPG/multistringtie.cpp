#include <vector> 
#ifndef NOTHREADS
#include "GThreads.h"
#include "GArgs.h"
#include "GStr.h"
#include "gff.h"
#include "GSam.h"
#include "GBitVec.h"
#include "tmerge.h"
#include "multist.h"
#include "unispg.h"
#include "infer_tranx.h"
// #include "unispg.h"

#include "time.h"
#include "GHashMap.hh"
#include "mode.hpp"
#endif

// Let's test single thread first!! 
#define NOTHREADS

//#define GMEMTRACE 1 
#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define VERSION "0.1.0"
#include <iostream>

#define USAGE "multiStringTie v" VERSION " usage:\n\n\
multiStringTie [CREATE_UNISPG / APPLY_UNISPG] <in.DOT ..> <in.bam ..> \
 [--mix] [--conservative] [--rf] [--fr]\n\
Assemble RNA-Seq alignments into potential transcripts with universal splice graph.\n\
Options:\n\
* Transcript CREATE_UNISPG mode: \n\
   multistringtie CREATE_UNISPG --multi [Options] { bam_list | sample1.bam ...}\n\
    --version : print just the version at stdout and exit\n\
 	-o output path/file name for the assembled transcripts GTF (default: stdout)\n\n\
* APPLY_UNISPG mode: \n\
   multistringtie APPLY_UNISPG [Options] { bam_list | sample1.bam ...}\n\
"

void processApplyOptions(GArgs& args);
char* sprintTime();
void unispg_readline(UnispgGp* &drec, bool &new_uni_spg_gp, bool &more_graph, int pre_refstart, int pre_refend);
void processBundleUnispg(BundleData* bundle, UnispgGp* unispg);

bool NoMoreBundles=false;
bool moreBundles(); //thread-safe retrieves NoMoreBundles
void noMoreBundles(); //sets NoMoreBundles to true

multiStringTieMode mode;
bool debugMode=false;
bool verbose=false;
GStr unispgdotfname_root; 
GStr unispgdotfname_pos; 
GStr unispgdotfname_neg; 


GStr outfname;
GStr out_dir;
GStr tmp_path;
GStr cram_ref; //reference genome FASTA for CRAM input
GStr tmpfname;
GStr genefname;
GStr traindir; // training directory for CDS option
bool keepTempFiles;
bool fr_strand=false;
bool rf_strand=false;
GStrSet<> excludeGseqs; //hash of chromosomes/contigs to exclude (e.g. chrM)
GffNames* gseqNames=NULL; //used as a dictionary for reference sequence names and ids
bool guided=false;
int refseqCount=0; // number of reference sequences found in the guides file
uint runoffdist=200;
double Num_Fragments=0; //global fragment counter (aligned pairs)
double Frag_Len=0;
double Cov_Sum=0;
int GeneNo=0; //-- global "gene" counter
float readthr=1;     // read coverage per bundle bp to accept it; // paper uses 3
float tpm_thr=1;
float fpkm_thr=1;
bool isunitig=true;
bool nomulti=false;
uint junctionsupport=10; // anchor length for junction to be considered well supported <- consider shorter??
uint sserror=25; // window arround splice sites that we use to generate consensus in case of long read data
int junctionthr=1; // number of reads needed to support a particular junction
float mcov=1; // fraction of bundle allowed to be covered by multi-hit reads paper uses 1
int mintranscriptlen=200; // minimum length for a transcript to be printed

TInputFiles bamreader;
DOTInputFile dotreader;
DOTInputFile dotreader_pos;
DOTInputFile dotreader_neg;

int sample_num = 1;

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


	} else if (mode == APPLY_UNISPG) {
		GArgs args(argc, argv,
		"debug;help;version;dot=;d=;bam;"
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


		// bool dot_is_open_pos = dotreader_pos.start(unispgdotfname_pos); //setup and open DOT input file
		// fprintf(stderr, "unispgdotfname_pos: %d\n", dot_is_open_pos);

		// bool dot_is_open_neg = dotreader_neg.start(unispgdotfname_neg); //setup and open DOT input file
		// fprintf(stderr, "unispgdotfname_neg: %d\n", dot_is_open_neg);

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








		/*****************************
		 * Processing Graphs & alignments: main algorithm
		 *****************************/
		/*****************************
		 * Iterating graphs
		 *****************************/
// 		while (more_graph_pos || more_graph_neg) {
// 			// After this, the current drec_pos is the newest positive graph.
// 			if (next_graph_pos) {
// 				if ((drec_pos=dotreader_pos.next())!=NULL) {
// 					fprintf(stderr, "***** drec_pos: %d - %d\n", drec_pos->get_refstart(), drec_pos->get_refend());
// 				} else {
// 					more_graph_pos = false;
// 					next_graph_pos = false;
// 				}
// 			}

// 			next_graph_neg = false;
// 			// After this, the current drec_pos is the newest negative graph.
// 			if (next_graph_neg) {
// 				if ((drec_neg=dotreader_neg.next())!=NULL) {
// 					fprintf(stderr, "***** drec_neg: %d - %d\n", drec_neg->get_refstart(), drec_neg->get_refend());
// 				} else {
// 					more_graph_neg = false;
// 					next_graph_neg = false;
// 				}
// 			}


// 			// /*****************************
// 			//  * Start processing the read.
// 			//  *****************************/
// 			// fprintf(stderr, "\t>>>>> Start processing read.\n");
// 			// fprintf(stderr, "\t>>>>> Read: %d - %d\n", r_start, r_end);

// 			// bool read_in_unispg_pos = true;
// 			// bool read_in_unispg_neg = true;
 
// 			// /*****************************
// 			//  * The read is on the positive strand 
// 			//  *****************************/
// 			// if ( xstrand=='+' ) {
// 			// 	if (r_end < drec_pos->get_refstart()) {
// 			// 		// it's out of the boundary. => go to the next read.
// 			// 		fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = false;
// 			// 	} else if (r_start < drec_pos->get_refstart() && r_end == drec_pos->get_refstart()) {
// 			// 		fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start < drec_pos->get_refstart() && r_end > drec_pos->get_refstart() && r_end < drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start < drec_pos->get_refstart() && r_end == drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start == drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end < drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start == drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end == drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start == drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end < drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end == drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 	} else if (r_start > drec_pos->get_refstart() && r_start == drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
// 			// 		next_graph_pos = false;
// 			// 		read_in_unispg_pos = true;
// 			// 		// int overlap_current = 0;
// 			// 		// overlap_current = drec->get_refend() - brec->start + 1;

// 			// 		// int overlap_next = 0;
// 			// 		// overlap_next = brec->end - drec->get_refend() + 1;

// 			// 		// if (overlap_current > overlap_next) {
// 			// 		// 	// The read belongs to the current processing bundle.
// 			// 		// } else {
// 			// 		// 	// The read belongs to the next bundle or not belongs to any bundles.
// 			// 		// 	read_in_unispg = false;
// 			// 		// }
// 			// 	} else if ( r_start > drec_pos->get_refend() ) {
// 			// 		fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
// 			// 		next_graph_pos = true;
// 			// 		read_in_unispg_pos = false;
// 			// 	} else {
// 			// 		GError("\tThis is a weird situation.!!! Check!!!\n");
// 			// 	}
// 			// }

// 			// /*****************************
// 			//  * The read is on the negative strand & it's out of the boundary.
// 			//  *    go to the next negative graph.
// 			//  *****************************/
// 			// if ( xstrand=='-' ) {
// 			// 	// if ( r_start > drec_neg->get_refend() ) {
// 			// 	// 	next_graph_neg = true;
// 			// 	// 	read_in_unispg_neg = false;
// 			// 	// }
// 			// }


// 	// 		/*****************************
// 	// 		 * Process the new bundle
// 	// 		 *****************************/
// 	// 		if (new_bundle) {
// 	// 			fprintf(stderr, "This is a new bundle\n");
// 	// 			fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
// 	// 			if (bundle->readlist.Count()>0) { // process reads in previous bundle
// 	// 				hashread.Clear();
// 	// 				// (readthr, junctionthr, mintranscriptlen are globals)
// 	// 				bundle->getReady(drec->get_refstart(), drec->get_refend());
// 	// 	#ifndef NOTHREADS
// 	// 				//push this in the bundle queue where it'll be picked up by the threads
// 	// 				DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
// 	// 				int qCount=0;
// 	// 				queueMutex.lock();
// 	// 				bundleQueue.Push(bundle);
// 	// 				bundleWork |= 0x02; //set bit 1
// 	// 				qCount=bundleQueue.Count();
// 	// 				queueMutex.unlock();
// 	// 				DBGPRINT2("##> bundleQueue.Count()=%d)\n", qCount);
// 	// 				//wait for a thread to pop this bundle from the queue
// 	// 				waitMutex.lock();
// 	// 				DBGPRINT("##> waiting for a thread to become available..\n");
// 	// 				while (threadsWaiting==0) {
// 	// 					haveThreads.wait(waitMutex);
// 	// 				}
// 	// 				waitMutex.unlock();
// 	// 				haveBundles.notify_one();
// 	// 				DBGPRINT("##> waitMutex unlocked, haveBundles notified, current thread yielding\n");
// 	// 				current_thread::yield();
// 	// 				queueMutex.lock();
// 	// 				DBGPRINT("##> queueMutex locked until bundleQueue.Count()==qCount\n");
// 	// 				while (bundleQueue.Count()==qCount) {
// 	// 					queueMutex.unlock();
// 	// 					DBGPRINT2("##> queueMutex unlocked as bundleQueue.Count()==%d\n", qCount);
// 	// 					haveBundles.notify_one();
// 	// 					current_thread::yield();
// 	// 					queueMutex.lock();
// 	// 					DBGPRINT("##> queueMutex locked again within while loop\n");
// 	// 				}
// 	// 				queueMutex.unlock();
// 	// 	#else //no threads
// 	// 				//Num_Fragments+=bundle->num_fragments;
// 	// 				//Frag_Len+=bundle->frag_len;
// 	// 				processBundleUnispg(bundle, drec);
// 	// 	#endif
// 	// 				// ncluster++; used it for debug purposes only
// 	// 			} //have alignments to process
// 	// 			else { //no read alignments in this bundle?
// 	// 	#ifndef NOTHREADS
// 	// 				dataMutex.lock();
// 	// 				DBGPRINT2("##> dataMutex locked for bundle #%d clearing..\n", bundle->idx);
// 	// 	#endif
// 	// 				bundle->Clear();
// 	// 	#ifndef NOTHREADS
// 	// 				dataClear.Push(bundle->idx);
// 	// 				DBGPRINT2("##> dataMutex unlocking as dataClear got pushed idx #%d\n", bundle->idx);
// 	// 				dataMutex.unlock();
// 	// 	#endif
// 	// 			} //nothing to do with this bundle
// 	// 			// if (chr_changed) {
// 	// 			// 	lastref=refseqName;
// 	// 			// 	lastref_id=gseq_id;
// 	// 			// 	// currentend=0;
// 	// 			// }
// 	// 			if (!more_alns) {
// 	// 				if (verbose) {
// 	// 	#ifndef NOTHREADS
// 	// 					GLockGuard<GFastMutex> lock(logMutex);
// 	// 	#endif
// 	// 					if (Num_Fragments) {
// 	// 						// printTime(stderr);
// 	// 						GMessage(" %g aligned fragments found.\n", Num_Fragments);
// 	// 					}
// 	// 					//GMessage(" Done reading alignments.\n");
// 	// 				}
// 	// 				noMoreBundles();
// 	// 				break;
// 	// 			}
// 	// #ifndef NOTHREADS
// 	// 			int new_bidx=waitForData(bundles);
// 	// 			if (new_bidx<0) {
// 	// 				//should never happen!
// 	// 				GError("Error: waitForData() returned invalid bundle index(%d)!\n",new_bidx);
// 	// 				break;
// 	// 			}
// 	// 			bundle=&(bundles[new_bidx]);
// 	// #endif
// 	// 			// currentstart=brec->start;
// 	// 			// currentend=brec->end;
// 	// 			bundle->refseq=lastref;
// 	// 			bundle->start=brec->start;
// 	// 			bundle->end= brec->end;
// 	// 			break;
// 	// 		}





// 			/*****************************
// 			 * Actually start processing reads.
// 			 *****************************/
// 			// if (brec != NULL) {
// 			// 	nh=brec->tag_int("NH");
// 			// 	if (nh==0) nh=1;
// 			// 	hi=brec->tag_int("HI");
// 			// 	// Old conditions to stop adding reads into a bundle
// 			// 	// if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
// 			// 	// 	new_bundle=true;
// 			// 	// }
// 			// 	// fprintf(stderr, "brec: %d\n", brec->mapped_len);
// 			// 	// if (brec->start > drec->get_refend()) {
// 			// 	// 	// fprintf(stderr, "brec: %d - %d\n", brec->start, brec->end);
// 			// 	// 	new_bundle=true;
// 			// 	// }
// 			// 	GReadAlnData alndata(brec, 0, nh, hi, tinfo);
// 			// 	if (xstrand=='+') alndata.strand=1;
// 			// 	else if (xstrand=='-') alndata.strand=-1;
// 			// 	//GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
// 			// 	//countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
// 			// 	//fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);

// 			// 	processRead(brec->start, brec->end, *bundle, hashread, alndata);
// 			// }
			


// 			/*****************************
// 			 * Iterating reads
// 			 *****************************/
// 			// while ( more_alns && ((read_in_unispg_pos && (xstrand=='+' || xstrand=='.') ) || (read_in_unispg_neg && (xstrand=='-' || xstrand=='.'))) ) {
// 			while (more_alns) {
// 				if ((brec=bamreader.next())!=NULL) {
// 					if (brec->isUnmapped()) continue;
// 					if (brec->start<1 || brec->mapped_len<10) {
// 						if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
// 								brec->name(), brec->start, brec->mapped_len);
// 						continue;
// 					}
// 		#ifdef DBG_ALN_DATA
// 					dbg_waln(brec);
// 		#endif
// 					r_start = brec->start; //start<end always!
// 					r_end = brec->end;

// 					refseqName=brec->refName();
// 					xstrand=brec->spliceStrand(); // tagged strand gets priority

// 					// fprintf(stderr, "** refseqName：%s\n", refseqName);

// 					/*****************************
// 					 * set strand if stranded library
// 					 *****************************/
// 					if(xstrand=='.' && (fr_strand || rf_strand)) {
// 						if(brec->isPaired()) { // read is paired
// 							if(brec->pairOrder()==1) { // first read in pair
// 								if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
// 								else xstrand='-';
// 							}
// 							else {
// 								if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
// 								else xstrand='+';
// 							}
// 						}
// 						else {
// 							if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
// 							else xstrand='-';
// 						}
// 					}





// 					/*****************************
// 					 * Start processing the read.
// 					 *****************************/
// 					fprintf(stderr, "\t>>>>> Start processing read.\n");
// 					fprintf(stderr, "\t>>>>> Read: %d - %d\n", r_start, r_end);

// 					bool read_in_unispg_pos = true;
// 					bool read_in_unispg_neg = true;
		
// 					/*****************************
// 					 * The read is on the positive strand 
// 					 *****************************/
// 					if ( xstrand=='+' ) {
// 						if (r_end < drec_pos->get_refstart()) {
// 							// it's out of the boundary. => go to the next read.
// 							fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = false;
// 						} else if (r_start < drec_pos->get_refstart() && r_end == drec_pos->get_refstart()) {
// 							fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start < drec_pos->get_refstart() && r_end > drec_pos->get_refstart() && r_end < drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start < drec_pos->get_refstart() && r_end == drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start == drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end < drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start == drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end == drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start == drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end < drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end == drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start > drec_pos->get_refstart() && r_start < drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 						} else if (r_start > drec_pos->get_refstart() && r_start == drec_pos->get_refend() && r_end > drec_pos->get_refend()) {
// 							fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
// 							next_graph_pos = false;
// 							read_in_unispg_pos = true;
// 							// int overlap_current = 0;
// 							// overlap_current = drec->get_refend() - brec->start + 1;

// 							// int overlap_next = 0;
// 							// overlap_next = brec->end - drec->get_refend() + 1;

// 							// if (overlap_current > overlap_next) {
// 							// 	// The read belongs to the current processing bundle.
// 							// } else {
// 							// 	// The read belongs to the next bundle or not belongs to any bundles.
// 							// 	read_in_unispg = false;
// 							// }
// 						} else if ( r_start > drec_pos->get_refend() ) {
// 							fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
// 							next_graph_pos = true;
// 							read_in_unispg_pos = false;
// 						} else {
// 							GError("\tThis is a weird situation.!!! Check!!!\n");
// 						}
// 					}

// 					/*****************************
// 					 * The read is on the negative strand & it's out of the boundary.
// 					 *    go to the next negative graph.
// 					 *****************************/
// 					if ( xstrand=='-' ) {
// 						// if ( r_start > drec_neg->get_refend() ) {
// 						// 	next_graph_neg = true;
// 						// 	read_in_unispg_neg = false;
// 						// }
// 					}

// 					total_read += 1;
// 					if (xstrand=='+') {
// 						pos_strand += 1;
// 					} else if (xstrand=='-') {
// 						neg_strand += 1;
// 					} else if (xstrand == '.') {
// 						unstrand += 1;
// 					}






// 	/*****************************
// 	 * Process the new bundle
// 	 *****************************/
// 	if (new_bundle) {
// 		fprintf(stderr, "This is a new bundle\n");
// 		fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
// 		if (bundle->readlist.Count()>0) { // process reads in previous bundle
// 			hashread.Clear();
// 			// (readthr, junctionthr, mintranscriptlen are globals)
// 			bundle->getReady(drec_pos->get_refstart(), drec_pos->get_refend());
// #ifndef NOTHREADS
// 			//push this in the bundle queue where it'll be picked up by the threads
// 			DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
// 			int qCount=0;
// 			queueMutex.lock();
// 			bundleQueue.Push(bundle);
// 			bundleWork |= 0x02; //set bit 1
// 			qCount=bundleQueue.Count();
// 			queueMutex.unlock();
// 			DBGPRINT2("##> bundleQueue.Count()=%d)\n", qCount);
// 			//wait for a thread to pop this bundle from the queue
// 			waitMutex.lock();
// 			DBGPRINT("##> waiting for a thread to become available..\n");
// 			while (threadsWaiting==0) {
// 				haveThreads.wait(waitMutex);
// 			}
// 			waitMutex.unlock();
// 			haveBundles.notify_one();
// 			DBGPRINT("##> waitMutex unlocked, haveBundles notified, current thread yielding\n");
// 			current_thread::yield();
// 			queueMutex.lock();
// 			DBGPRINT("##> queueMutex locked until bundleQueue.Count()==qCount\n");
// 			while (bundleQueue.Count()==qCount) {
// 				queueMutex.unlock();
// 				DBGPRINT2("##> queueMutex unlocked as bundleQueue.Count()==%d\n", qCount);
// 				haveBundles.notify_one();
// 				current_thread::yield();
// 				queueMutex.lock();
// 				DBGPRINT("##> queueMutex locked again within while loop\n");
// 			}
// 			queueMutex.unlock();
// #else //no threads
// 			//Num_Fragments+=bundle->num_fragments;
// 			//Frag_Len+=bundle->frag_len;
// 			processBundleUnispg(bundle, drec_pos);
// #endif
// 			// ncluster++; used it for debug purposes only
// 		} //have alignments to process
// 		else { //no read alignments in this bundle?
// #ifndef NOTHREADS
// 			dataMutex.lock();
// 			DBGPRINT2("##> dataMutex locked for bundle #%d clearing..\n", bundle->idx);
// #endif
// 			bundle->Clear();
// #ifndef NOTHREADS
// 			dataClear.Push(bundle->idx);
// 			DBGPRINT2("##> dataMutex unlocking as dataClear got pushed idx #%d\n", bundle->idx);
// 			dataMutex.unlock();
// #endif
// 		} //nothing to do with this bundle
// 		// if (chr_changed) {
// 		// 	lastref=refseqName;
// 		// 	lastref_id=gseq_id;
// 		// 	// currentend=0;
// 		// }
// 		if (!more_alns) {
// 			if (verbose) {
// #ifndef NOTHREADS
// 				GLockGuard<GFastMutex> lock(logMutex);
// #endif
// 				if (Num_Fragments) {
// 					// printTime(stderr);
// 					GMessage(" %g aligned fragments found.\n", Num_Fragments);
// 				}
// 				//GMessage(" Done reading alignments.\n");
// 			}
// 			noMoreBundles();
// 			break;
// 		}
// #ifndef NOTHREADS
// 		int new_bidx=waitForData(bundles);
// 		if (new_bidx<0) {
// 			//should never happen!
// 			GError("Error: waitForData() returned invalid bundle index(%d)!\n",new_bidx);
// 			break;
// 		}
// 		bundle=&(bundles[new_bidx]);
// #endif
// 		// currentstart=brec->start;
// 		// currentend=brec->end;
// 		bundle->refseq=lastref;
// 		bundle->start=brec->start;
// 		bundle->end= brec->end;
// 		break;
// 	}








// 					nh=brec->tag_int("NH");
// 					if (nh==0) nh=1;
// 					hi=brec->tag_int("HI");

// 					// Old conditions to stop adding reads into a bundle
// 					// if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
// 					// 	new_bundle=true;
// 					// }
// 					// fprintf(stderr, "brec: %d\n", brec->mapped_len);
// 					// if (brec->start > drec->get_refend()) {
// 					// 	// fprintf(stderr, "brec: %d - %d\n", brec->start, brec->end);
// 					// 	new_bundle=true;
// 					// }


// 					GReadAlnData alndata(brec, 0, nh, hi, tinfo);
// 					if (xstrand=='+') alndata.strand=1;
// 					else if (xstrand=='-') alndata.strand=-1;
// 					//GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
// 					//countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
// 					//fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);

// 					processRead(brec->start, brec->end, *bundle, hashread, alndata);





// 				} else {
// 					more_alns=false;
// 				}
// 			}

// 			fprintf(stderr, "total read num: %d\n", total_read);
// 			fprintf(stderr, "positve read num: %d\n", pos_strand);
// 			fprintf(stderr, "negative read num: %d\n", neg_strand);
// 			fprintf(stderr, "unstranded read num: %d\n", unstrand);

// 		}


		while (more_graph) {
			if ((drec=dotreader.next())!=NULL) {

				while (more_alns) {
					if ((brec=bamreader.next())!=NULL) {

						prev_brec= brec;
					} else {
						more_alns=false;
						new_bundle=true; //fake a new start (end of last bundle)
					}

					//  If reads overlap with the previous unispg
					//    skip if it does not overlap


					// Process the previous unispg bundle.
				}
				prev_drec = drec;
			} else {
				more_graph=false;
			}
		}










		while (more_graph) {
			// unispg_readline(drec, new_uni_spg_gp, more_graph, pre_refstart, pre_refend) {}

			if ((drec=dotreader.next())!=NULL) {
				fprintf(stderr, "***** drec: %d - %d\n", drec->get_refstart(), drec->get_refend());

				// drec->get_no2gnode();
				// GVec<bool>* is_passed_source = new GVec<bool>(sample_num-1, false);
				// GVec<float>* cov_source = new GVec<float>(sample_num-1, 0.0f);
				// GVec<float>* capacity_source = new GVec<float>(sample_num-1, 0.0f);
				// CGraphnodeUnispg* source = new CGraphnodeUnispg(sample_num, 0, 0, new_unispg_nodeid, is_passed_source, cov_source, capacity_source, true, 0, 0, 0, true, 0, 0);


				// Move on to the next graph.
				if (r_start > drec->get_refend()) {
					fprintf(stderr, "** Prev read start is bigger than graph coordinate\n");
					continue;
				}

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

						/*****************************
						 * set strand if stranded library
						 *****************************/
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


						/*****************************
						 * Found the overlapping between a read & the unispg
						 *****************************/
						if (xstrand == '+') {
							/*****************************
							 ** Step 1-2: Check whether reads are in the range of the graph.
							*****************************/
							// fprintf(stderr, "Process read (pre: %d - %d ;  now: %d - %d)!!!\n", pre_refstart, pre_refend, brec->start, brec->end);
							fprintf(stderr, "\t>>>>> brec: %d - %d\n", brec->start, brec->end);
							if (brec->end <  drec->get_refstart()) {
								// ----------   |(s).................(e)|
								// The read is outside the current bundle => skipped!
								fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
								process_read = false;
								continue;
								// read_in_unispg = false;
								// new_bundle = true;
							}
							if (brec->start < drec->get_refstart() && brec->end == drec->get_refstart()) {
								// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
								fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
							} else if (brec->start < drec->get_refstart() && brec->end > drec->get_refstart()) {
								// ----------|(s).................(e)|   or   -----|(s)-----............(e)|
								fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
							} else if (brec->start < drec->get_refstart() && brec->end == drec->get_refstart()) {
								// -----|(s)---------(e)|
								fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
							} else if (brec->start == drec->get_refstart() && brec->start < drec->get_refend() && brec->end < drec->get_refend()) {
								// |(s)----------.................(e)|   or   |(s)....----------........(e)|
								fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
							} else if (brec->start == drec->get_refstart() && brec->start < drec->get_refend() && brec->end == drec->get_refend()) {
								// |(s)----------.................(e)|   or   |(s)....----------........(e)|
								fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
							} else if (brec->start == drec->get_refstart() && brec->start < drec->get_refend() && brec->end > drec->get_refend()) {
								// |(s)----------.................(e)|   or   |(s)....----------........(e)|
								fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
							} else if (brec->start > drec->get_refstart() && brec->start < drec->get_refend() && brec->end < drec->get_refend()) {
								// |(s)----------.................(e)|   or   |(s)....----------........(e)|
								fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
							} else if (brec->start > drec->get_refstart() && brec->start < drec->get_refend() && brec->end == drec->get_refend()) {
								// |(s)----------.................(e)|   or   |(s)....----------........(e)|
								fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
							} else if (brec->start > drec->get_refstart() && brec->start < drec->get_refend() && brec->end > drec->get_refend()) {
								// |(s)----------.................(e)|   or   |(s)....----------........(e)|
								fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
							} else if (brec->start > drec->get_refstart() && brec->start < drec->get_refend() && brec->end > drec->get_refend()) {
								// |(s)...............------(e)|-----    or   |(s).................(e)|----------   
								// The overlapping with the current processing bundle.
								fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
							} else if (brec->start > drec->get_refstart() && brec->start == drec->get_refend() && brec->end > drec->get_refend()) {
								// |(s)...............------(e)|-----    or   |(s).................(e)|----------   
								// The overlapping with the current processing bundle.
								fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
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
							} else if (brec->start > drec->get_refend()) {
								fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
								read_in_unispg = false;
								new_bundle = true;
								process_read = false;
							}
						// } else if (xstrand == '-' || xstrand == '.') {
						} else if (xstrand == '-') {
							continue;
						}








					} else { //no more alignments
						more_alns=false;
						new_bundle=true; //fake a new start (end of last bundle)
					}










					/*****************************
					 * Condition to start processing reads in the previous bundle
					 *****************************/
					// if (new_bundle || chr_changed) {
					if (new_bundle) {
						fprintf(stderr, "This is a new bundle\n");
						fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
						hashread.Clear();
				
				
				
						if (bundle->readlist.Count()>0) { // process reads in previous bundle
							// (readthr, junctionthr, mintranscriptlen are globals)
							bundle->getReady(drec->get_refstart(), drec->get_refend());

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
							processBundleUnispg(bundle, drec);
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
						bundle->end=brec->end;
						break;
					} //<---- new bundle started
				
				















					if (process_read && brec!=NULL) {
						/*****************************
						 * Actually start processing reads.
						 *****************************/
						// if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
						pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
						// chr_changed=(lastref.is_empty() || lastref!=refseqName);
						// if (chr_changed) {
						// 	skipGseq=excludeGseqs.hasKey(refseqName);
						// 	gseq_id=gseqNames->gseqs.addName(refseqName);
						// 	if (guided) {
						// 		if (gseq_id>=refseqCount) {
						// 			if (verbose)
						// 				GMessage("WARNING: no reference transcripts found for genomic sequence \"%s\"! (mismatched reference names?)\n",
						// 					refseqName);
						// 		}
						// 		else no_ref_used=false;
						// 	}
						// 	if (alncounts.Count()<=gseq_id) {
						// 		alncounts.Resize(gseq_id+1);
						// 	}
						// 	else if (alncounts[gseq_id]>0)
						// 			GError("%s\nAlignments (%d) already found for %s !\n",
						// 				ERR_BAM_SORT, alncounts[gseq_id], refseqName);
						// 	prev_pos=0;
						// }

						if (pos<prev_pos) GError("%s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
							ERR_BAM_SORT, brec->name(), brec->start,  pos, refseqName, prev_pos);
						prev_pos=pos;
						// if (skipGseq) continue;
						alncounts[gseq_id]++;
						nh=brec->tag_int("NH");
						if (nh==0) nh=1;
						hi=brec->tag_int("HI");

						// Old conditions to stop adding reads into a bundle
						// if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
						// 	new_bundle=true;
						// }
						// fprintf(stderr, "brec: %d\n", brec->mapped_len);
						// if (brec->start > drec->get_refend()) {
						// 	// fprintf(stderr, "brec: %d - %d\n", brec->start, brec->end);
						// 	new_bundle=true;
						// }


						GReadAlnData alndata(brec, 0, nh, hi, tinfo);
						// bool ovlpguide=bundle->evalReadAln(alndata, xstrand);

						/*****************************
						** Step 4: in eonly case consider read only if it overlaps guide
						**     check for overlaps with ref transcripts which may set xstrand
						*****************************/
						// eonly: for mergeMode includes estimated coverage sum in the merged transcripts
						// fprintf(stderr, "ovlpguide: %d\n", ovlpguide);
						// if(ovlpguide) {
						if (xstrand=='+') {
							alndata.strand=1;
							pos_strand += 1;
						} else if (xstrand=='-') {
							alndata.strand=-1;
							neg_strand += 1;
						} else if (xstrand=='.') {
							unstrand += 1;
						}
						//GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
						//countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
						//fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);

						if (xstrand=='+' || xstrand=='.') {
							fprintf(stderr, ">> xstrand: %c\n", xstrand);
							processRead(brec->start, brec->end, *bundle, hashread, alndata);
							processed_read += 1;
						}
						// }

						// if (new_bundle && bundle->readlist.Count() == 1) {
						// 	break;
						// }	
					}















					prev_brec= brec;
				} //for each read alignment
				prev_drec = drec;
			} else {
				// fprintf(stderr, "No more dot graph!!!\n");
				more_graph=false;
				// new_uni_spg_gp = true; //fake a new start (end of the last universal splice graph.)
			}
		}
		fprintf(stderr, "total read num: %d\n", total_read);
		fprintf(stderr, "processed read num: %d\n", processed_read);
		fprintf(stderr, "positve read num: %d\n", pos_strand);
		fprintf(stderr, "negative read num: %d\n", neg_strand);
		fprintf(stderr, "unstranded read num: %d\n", unstrand);
		delete brec;
		bamreader.stop(); //close all BAM files
		delete drec;
		dotreader.stop(); //close all DOT files
		// dotreader_pos.stop(); //close all DOT files
		// dotreader_neg.stop(); //close all DOT files
	}
    return 0;
}





void processApplyOptions(GArgs& args) {

   fprintf(stderr, "Inside 'processApplyOptions'\n");
	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n",VERSION);
	   exit(0);
	}	

	tmpfname=args.getOpt('o');
	GStr s=args.getOpt("dot");
	// if (s.is_empty()) {
	//    s=args.getOpt("d");
	//    fprintf(stderr, "It's empty%s \n", s.chars());
	// }
	if (!s.is_empty()) {
		unispgdotfname_root=s;
	}

	if (args.getOpt("fr")) {
		fr_strand=true;
	}
	if (args.getOpt("rf")) {
		rf_strand=true;
		if(fr_strand) GError("Error: --fr and --rf options are incompatible.\n");
	}


   fprintf(stderr, "%s \n", unispgdotfname_root.chars());




	const char* ifn=NULL;
	while ( (ifn=args.nextNonOpt())!=NULL) {
		//input alignment files
		bamreader.Add(ifn);
	}
	//deferred creation of output path
	outfname="stdout";
	out_dir="./";
	if (!tmpfname.is_empty() && tmpfname!="-") {
		if (tmpfname[0]=='.' && tmpfname[1]=='/')
			tmpfname.cut(0,2);
		outfname=tmpfname;
		int pidx=outfname.rindex('/');
		if (pidx>=0) {//path given
			out_dir=outfname.substr(0,pidx+1);
			tmpfname=outfname.substr(pidx+1);
		}
	} else { // stdout
		tmpfname=outfname;
		char *stime=sprintTime();
		tmpfname.tr(":","-");
		tmpfname+='.';
		tmpfname+=stime;
	}
	if (out_dir!="./") {
		if (fileExists(out_dir.chars())==0) {
			//directory does not exist, create it
			if (Gmkdir(out_dir.chars()) && !fileExists(out_dir.chars())) {
				GError("Error: cannot create directory %s!\n", out_dir.chars());
			}
		}
	}

	if (!genefname.is_empty()) {
		if (genefname[0]=='.' && genefname[1]=='/')
			genefname.cut(0,2);
		//attempt to create the gene abundance path
		GStr genefdir("./");
		int pidx=genefname.rindex('/');
		if (pidx>=0) { //get the path part
			genefdir=genefname.substr(0,pidx+1);
			//genefname=genefname.substr(pidx+1);
		}
		if (genefdir!="./") {
			if (fileExists(genefdir.chars())==0) {
				//directory does not exist, create it
				if (Gmkdir(genefdir.chars()) || !fileExists(genefdir.chars())) {
					GError("Error: cannot create directory %s!\n", genefdir.chars());
				}
			}
		}

	}

	{ //prepare temp path
		GStr stempl(out_dir);
		stempl.chomp('/');
		stempl+="/tmp_XXXXXX";
		char* ctempl=Gstrdup(stempl.chars());
		Gmktempdir(ctempl);
		tmp_path=ctempl;
		tmp_path+='/';
		GFREE(ctempl);
	}

	tmpfname=tmp_path+tmpfname;
   	fprintf(stderr, "%s \n", tmpfname.chars());
}

//----------------------------------------
char* sprintTime() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}

// void unispg_readline(UnispgGp* &drec, bool &new_uni_spg_gp, bool &more_graph, int pre_refstart, int pre_refend) {
// 	if ((drec=dotreader.next())!=NULL) {
// 		if (pre_refstart != drec->get_refstart() && pre_refend != drec->get_refend()) {
// 			/*****************************************
// 			 ** It is a new universal splice graph group.
// 			 *****************************************/
// 			new_uni_spg_gp = true;
// 		} 
// 		// fprintf(stderr, "Keep reading %d (pre: %d - %d ;  now: %d - %d)!!!\n", new_uni_spg_gp, pre_refstart, pre_refend, drec->get_refstart(), drec->get_refend());


// 		// int uni_refstart = drec->get_refstart();
// 		// int uni_refend = drec->get_refend();
// 		// int uni_s = drec->get_s();
// 		// int uni_g = drec->get_g_idx();
// 		// int uni_graphno = drec->get_graphno();
// 		// int uni_edgeno = drec->get_edgeno();
// 		// GPVec<CGraphnode>* uni_no2gnode = drec->get_no2gnode(); 

// 		// for (int i = 1; i < uni_graphno-1; i++) {
// 		// 	fprintf(stderr, "uni_no2gnode[s][g][i]: %d  start: %d   end: %d \n", uni_no2gnode[0][i]->nodeid, uni_no2gnode[0][i]->start, uni_no2gnode[0][i]->end);
// 		// }


// 	} else {
// 		// fprintf(stderr, "No more dot graph!!!\n");
// 		more_graph=false;
// 		new_uni_spg_gp = true; //fake a new start (end of the last universal splice graph.)
// 	}
// }

void processBundleUnispg(BundleData* bundle, UnispgGp* unispg) {
	fprintf(stderr, "Inside processBundleUnispg!\n");
	if (verbose) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(logMutex);
#endif
		// printTime(stderr);
		GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d junctions, %d guides] begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads, bundle->readlist.Count(), bundle->junction.Count(),
                bundle->keepguides.Count());
#ifdef GMEMTRACE
			double vm,rsm;
			get_mem_usage(vm, rsm);
			GMessage("\t\tstart memory usage: %6.1fMB\n",rsm/1024);
			if (rsm>maxMemRS) {
				maxMemRS=rsm;
				maxMemVM=vm;
				maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
			}
#endif
	}
#ifdef B_DEBUG
	for (int i=0;i<bundle->keepguides.Count();++i) {
		GffObj& t=*(bundle->keepguides[i]);
		RC_TData* tdata=(RC_TData*)(t.uptr);
		fprintf(dbg_out, ">%s (t_id=%d) %s%c %d %d\n", t.getID(), tdata->t_id, t.getGSeqName(), t.strand, t.start, t.end );
		for (int fe=0;fe < tdata->t_exons.Count(); ++fe) {
			RC_Feature& exoninfo = *(tdata->t_exons[fe]);
			fprintf(dbg_out, "%d\texon\t%d\t%d\t%c\t%d\t%d\n", exoninfo.id, exoninfo.l, exoninfo.r,
					    exoninfo.strand, exoninfo.rcount, exoninfo.ucount);
			if (! (exoninfo==*(bundle->rc_data->guides_RC_exons->Get(exoninfo.id-1))))
				 GError("exoninfo with id (%d) not matching!\n", exoninfo.id);
		}
		for (int fi=0;fi < tdata->t_introns.Count(); ++fi) {
			RC_Feature& introninfo = *(tdata->t_introns[fi]);
			fprintf(dbg_out, "%d\tintron\t%d\t%d\t%c\t%d\t%d\n", introninfo.id, introninfo.l, introninfo.r,
					introninfo.strand, introninfo.rcount, introninfo.ucount);
			if (! (introninfo==*(bundle->rc_data->guides_RC_introns->Get(introninfo.id-1))))
				 GError("introninfo with id (%d) not matching!\n", introninfo.id);
		}
		//check that IDs are properly assigned
		if (tdata->t_id!=(uint)t.udata) GError("tdata->t_id(%d) not matching t.udata(%d)!\n",tdata->t_id, t.udata);
		if (tdata->t_id!=bundle->rc_data->guides_RC_tdata->Get(tdata->t_id-1)->t_id)
			 GError("tdata->t_id(%d) not matching rc_data[t_id-1]->t_id (%d)\n", tdata->t_id, bundle->rc_data->g_tdata[tdata->t_id-1]->t_id);

	}
#endif

	infer_transcripts_unispg(bundle, unispg);

	if (bundle->pred.Count()>0) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		// GeneNo=printResults(bundle, GeneNo, bundle->refseq);
	}

	if (bundle->num_fragments) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(countMutex);
		#endif
		Num_Fragments+=bundle->num_fragments;
		Frag_Len+=bundle->frag_len;
		Cov_Sum+=bundle->sum_cov;
	}

	if (verbose) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(logMutex);
#endif
	  /*
	  SumReads+=bundle->sumreads;
	  SumFrag+=bundle->sumfrag;
	  NumCov+=bundle->num_cov;
	  NumReads+=bundle->num_reads;
	  NumFrag+=bundle->num_frag;
	  NumFrag3+=bundle->num_fragments3;
	  SumFrag3+=bundle->sum_fragments3;
	  fprintf(stderr,"Number of fragments in bundle: %g with length %g\n",bundle->num_fragments,bundle->frag_len);
	  fprintf(stderr,"Number of fragments in bundle: %g with sum %g\n",bundle->num_fragments,bundle->frag_len);
	  */
		// printTime(stderr);
		GMessage("^bundle %s:%d-%d done (%d processed potential transcripts).\n",bundle->refseq.chars(),
				bundle->start, bundle->end, bundle->pred.Count());
#ifdef GMEMTRACE
		double vm,rsm;
		get_mem_usage(vm, rsm);
		GMessage("\t\tfinal memory usage: %6.1fMB\n",rsm/1024);
		if (rsm>maxMemRS) {
			maxMemRS=rsm;
			maxMemVM=vm;
			maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
		}
#endif
	}
	bundle->Clear();
}

void noMoreBundles() {
#ifndef NOTHREADS
		bamReadingMutex.lock();
		NoMoreBundles=true;
		bamReadingMutex.unlock();
		queueMutex.lock();
		bundleWork &= ~(int)0x01; //clear bit 0;
		queueMutex.unlock();
		bool areThreadsWaiting=true;
		do {
		  waitMutex.lock();
		   areThreadsWaiting=(threadsWaiting>0);
		  waitMutex.unlock();
		  if (areThreadsWaiting) {
		    DBGPRINT("##> NOTIFY ALL workers: no more data!\n");
		    haveBundles.notify_all();
		    current_thread::sleep_for(1);
		    waitMutex.lock();
		     areThreadsWaiting=(threadsWaiting>0);
		    waitMutex.unlock();
		    current_thread::sleep_for(1);
		  }
		} while (areThreadsWaiting); //paranoid check that all threads stopped waiting
#else
	  NoMoreBundles=true;
#endif
}