#pragma once
#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "rlink.h"
#include "t_record.h"
#include "dot_record.h"
#include "unispg.h"
#include "helper.h"
#include "definitions.h"

#include "APPLY_UNISPG/processOptions_A.h"
#include "CREATE_UNISPG/processOptions_C.h"

#include "APPLY_UNISPG/processBundle_A.h"
#include "CREATE_UNISPG/processBundle_C.h"

#include "APPLY_UNISPG/unispg_A.h"
#include "CREATE_UNISPG/unispg_C.h"

// Include fstream to quickly write out the ratio!!!
#include <iostream>
#include <fstream>
#include <vector> 
#include "time.h"

/*******************************************
 ** Argument parsing parameters.
 *******************************************/
bool debugMode=false; // "debug" or "D" tag.
bool verbose=false; // "verbose" / "v" tag.
bool ballgown=false; // "B" tag.
GStr ballgown_dir; // "b" tag.
bool viral=false; // "viral" tag.
bool mixedMode=false; // "mix" tag. both short and long read data alignments are provided
bool mergeMode = false; // "merge" tag. For running StringTie Merge.
bool multiMode=false; // "multi" tag.
bool graph_bed=false; // "graph_bed" tag.
bool keepTempFiles; // "keeptmp" tag.
GFastaDb* gfasta=NULL; // "rseq" or "S" tag.
GStr ptff; // "ptf" tag. (point features)
bool fr_strand=false; // "fr" tag.
bool rf_strand=false; // "rf" tag.
bool includesource=true; // "z" tag.
bool retained_intron=false; // "i" tag. set by parameter -i for merge option
bool trim=true; // "t" tag. 
bool eonly=false; // "e" tag. for mergeMode includes estimated coverage sum in the merged transcripts
bool nomulti=false; // "u" tag.
bool longreads=false; // "L" tag.
bool rawreads=false; // "R" tag.
GStrSet<> excludeGseqs; // "x" tag. hash of chromosomes/contigs to exclude (e.g. chrM)
bool guided=false; // "G" tag.
GStr guidegff; // "G" tag
GStr label("STRG"); // "l" tag.
int mintranscriptlen=200; // "m" tag. minimum length for a transcript to be printed
uint junctionsupport=10; // "a" tag. anchor length for junction to be considered well supported <- consider shorter??
int junctionthr=1; // "j" tag. number of reads needed to support a particular junction
float singlethr=4.75; // "s" tag. coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
float readthr=1; // "c" tag. read coverage per bundle bp to accept it; // paper uses 3
float isofrac=0.01; // "f" tag. 
int num_cpus=1; // "p" tag.
uint bundledist=50;  // "g" tag. reads at what distance should be considered part of separate bundles
uint runoffdist=200; // threshold for 'bundledist'
float mcov=1; // "M" tag. fraction of bundle allowed to be covered by multi-hit reads paper uses 1
bool geneabundance=false; // "A" tag.
uint sserror=25; // "E" tag. window arround splice sites that we use to generate consensus in case of long read data
float fpkm_thr=1; // "F" tag.
float tpm_thr=1; // "T" tag.
bool isunitig=true; // "U" tag.
bool enableNames=false; // "E" tag (removed)
FILE* c_out=NULL; // "C" tag.
FILE* f_out=NULL; // get from "outfname" param.


/*******************************************
 ** File-related parameters.
 *******************************************/
// Dot file => for reading in DOT format.
GStr unispgdotfname_root; 
GStr unispgdotfname_pos; 
GStr unispgdotfname_neg; 
// Annotation file => outputing in GFF format.
GStr outfname;
GStr outfname_prefix;
GStr out_dir;
GStr tmp_path;
GStr cram_ref; //"ref" / "cram-ref" tag. Reference genome FASTA for CRAM input
GStr tmpfname; // "o" tag.
GStr genefname;
GStr traindir; // "cds" tag. training directory for CDS option (removed)
// Ratio file => for coverage comparison & visualization.
ofstream cov_file_pos;
ofstream cov_file_neg;
ofstream cov_file_pos_norm;
ofstream cov_file_neg_norm;

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
int allowed_nodes=1000;

/*******************************************
 ** Reader parameters.
 *******************************************/
TInputFiles bamreader;
DOTInputFile dotreader;
DOTInputFile dotreader_pos;
DOTInputFile dotreader_neg;



/*******************************************
 ** CREATE_UNISPG specific parameters.
 *******************************************/
UnispgGp_CREATE* unispg_gp;
bool universal_splice_graph = false;
FILE* uinigraph_out = NULL;
GStr unigraphfname; 
GStr plot_dir;

/*****************************
 * Declaring DOT file.
 *****************************/
GVec<FILE*> pos_dot_vec = NULL;
GVec<GStr> pos_dotfname_vec; 
GVec<FILE*> neg_dot_vec = NULL;
GVec<GStr> neg_dotfname_vec; 
FILE* dot = NULL;
GStr dotfname; 
GVec<FILE*>* dot_vec[2];
GVec<GStr>* dotfname_vec[2]; 

/*****************************
 * Declaring BED file.
 *****************************/
GVec<FILE*>* node_lclg_bed_vec[2];
GVec<GStr>* nodelclgfname_vec[2]; 
GVec<FILE*>* edge_lclg_bed_vec[2];
GVec<GStr>* edgelclgfname_vec[2]; 

GVec<FILE*>* node_novp_bed_vec[2];
GVec<GStr>* nodenovpfname_vec[2]; 
GVec<FILE*>* edge_novp_bed_vec[2];
GVec<GStr>* edgenovpfname_vec[2]; 

GVec<FILE*>* node_unispg_bed_vec[2];
GVec<GStr>* nodeunispgfname_vec[2]; 
GVec<FILE*>* edge_unispg_bed_vec[2];
GVec<GStr>* edgeunispgfname_vec[2]; 

FILE* node_unispg_unstrand_bed; 
GStr nodeunispgfname_unstrand; 
// These are for temporary file holding
FILE* node_cov_bed = NULL;
GStr nodecovfname; 
FILE* edge_cov_bed = NULL;
GStr edgecovfname; 

/*****************************
 * Declaring reference related data structure
 *****************************/
GVec<GRefData> refguides; // plain vector with transcripts for each chromosome
GArray<GRefPtData> refpts(true, true); // sorted,unique array of refseq point-features data

/*****************************
 * Declaring Ballgown related data structure
 *   table indexes for Ballgown Raw Counts data (-B/-b option)
 *****************************/
GPVec<RC_TData> guides_RC_tdata(true); //raw count data or other info for all guide transcripts
GPVec<RC_Feature> guides_RC_exons(true); //raw count data for all guide exons
GPVec<RC_Feature> guides_RC_introns(true);//raw count data for all guide introns




/*******************************************
 ** Multithread parameters.
 *******************************************/
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
		/*******************************************
		 ** Process arguments.
		 *******************************************/
		GArgs args(argc, argv,
		"debug;help;version;ballgown;viral;conservative;mix;merge;multi;graph_bed;ref=;cram-ref=cds=;keeptmp;rseq=;ptf=;fr;rf;"
		"exclude=zihvteuLRx:n:j:s:D:G:C:S:l:m:o:a:c:f:p:g:M:Bb:A:E:F:T:");
		args.printError(USAGE, true);
		processCreateOptions(args);

		GVec<int> alncounts(30); //keep track of the number of read alignments per chromosome [gseq_id]

		int bamcount=bamreader.start(); //setup and open input files
		fprintf(stderr, "&& bamcount number: %d\n", bamcount);

		/*******************************************
		 ** Initiating all file-related parameters.
		 *******************************************/
		for(int s=0;s<2;s++) { // skip neutral bundles -> those shouldn't have junctions
			dot_vec[s] = new GVec<FILE*>[bamcount];
			dotfname_vec[s] = new GVec<GStr>[bamcount];
			if (graph_bed) {
				// Local graph (lclg)
				node_lclg_bed_vec[s] = new GVec<FILE*>[bamcount];
				nodelclgfname_vec[s] = new GVec<GStr>[bamcount];
				edge_lclg_bed_vec[s] = new GVec<FILE*>[bamcount];
				edgelclgfname_vec[s] = new GVec<GStr>[bamcount];
				// Nonoverlap graph (novp)
				node_novp_bed_vec[s] = new GVec<FILE*>[bamcount];
				nodenovpfname_vec[s] = new GVec<GStr>[bamcount];
				edge_novp_bed_vec[s] = new GVec<FILE*>[bamcount];
				edgenovpfname_vec[s] = new GVec<GStr>[bamcount];
				// Universal splice graph (unispg)
				node_unispg_bed_vec[s] = new GVec<FILE*>[bamcount];
				nodeunispgfname_vec[s] = new GVec<GStr>[bamcount];
				edge_unispg_bed_vec[s] = new GVec<FILE*>[bamcount];
				edgeunispgfname_vec[s] = new GVec<GStr>[bamcount];
			}
		}

	#ifndef GFF_DEBUG
		if (bamcount<1) {
			GError("%sError: no input files provided!\n",USAGE);
		}
	#endif
	#ifdef DEBUGPRINT
		verbose=true;
	#endif
		const char* ERR_BAM_SORT="\nError: the input alignment file is not sorted!\n";

		unispg_gp = new UnispgGp_CREATE();
		plot_dir = outfname.copy();		
		if (outfname.endsWith(".gtf")) {
			plot_dir.chomp(".gtf");
		}
		if (fileExists(plot_dir.chars())==0) {
			//directory does not exist, create it
			if (Gmkdir(plot_dir.chars()) && !fileExists(plot_dir.chars())) {
				GError("Error: cannot create directory %s!\n", plot_dir.chars());
			}
		}
		outfname_prefix = outfname.copy();		
		if (outfname.endsWith(".gtf")) {
			outfname_prefix.chomp(".gtf");
		}
		fprintf(stderr, "outfname: %s\n", outfname.chars());
		fprintf(stderr, "out_dir: %s\n", out_dir.chars());
		fprintf(stderr, "tmp_path: %s\n", tmp_path.chars());
		fprintf(stderr, "cram_ref: %s\n", cram_ref.chars());
		fprintf(stderr, "tmpfname: %s\n", tmpfname.chars());
		fprintf(stderr, "genefname: %s\n", genefname.chars());
		fprintf(stderr, "traindir: %s\n", traindir.chars());


		/*******************************************
		 ** read guiding transcripts from input gff file
		 *******************************************/
		if(guided) { 
			if (verbose) {
				printTime(stderr);
				GMessage(" Loading reference annotation (guides)..\n");
			}
			FILE* f=fopen(guidegff.chars(),"r");
			if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
				guidegff.chars());
			//                transcripts_only    sort by location?
			GffReader gffr(f,       true,             true); //loading only recognizable transcript features
			gffr.setRefAlphaSorted(); //alphabetical sorting of refseq IDs
			gffr.showWarnings(verbose);
			//        keepAttrs    mergeCloseExons   noExonAttrs
			gffr.readAll(false,          true,        true);
			//the list of GffObj is in gffr.gflst, sorted by chromosome and start-end coordinates
			//collect them in other data structures, if it's kept for later call gffobj->isUsed(true)
			// (otherwise it'll be deallocated when gffr is destroyed due to going out of scope)
			refseqCount=gffr.gseqtable.Count();
			if (refseqCount==0 || gffr.gflst.Count()==0) {
				GError("Error: could not any valid reference transcripts in %s (invalid GTF/GFF file?)\n",
						guidegff.chars());
			}
			refguides.setCount(refseqCount); //maximum gseqid
			uint c_tid=0;
			uint c_exon_id=0;
			uint c_intron_id=0;
			GList<RC_Feature> uexons(true, false, true); //sorted, free items, unique
			GList<RC_Feature> uintrons(true, false, true);
			//assign unique transcript IDs based on the sorted order
			int last_refid=-1;
			bool skipGseq=false;
			for (int i=0;i<gffr.gflst.Count();i++) {
				GffObj* m=gffr.gflst[i];
				if (last_refid!=m->gseq_id) {
					//chromosome switch
					if (ballgown) { //prepare memory storage/tables for all guides on this chromosome
						uexons.Clear();
						uintrons.Clear();
					}
					last_refid=m->gseq_id;
					skipGseq=excludeGseqs.hasKey(m->getGSeqName());
				}
				//sanity check: make sure there are no exonless "genes" or other
				if (skipGseq) continue;
				if (m->exons.Count()==0) {
					if (verbose)
						GMessage("Warning: exonless GFF %s feature with ID %s found, added implicit exon %d-%d.\n",
								m->getFeatureName(), m->getID(), m->start, m->end);
					m->addExon(m->start, m->end); //should never happen!
				}
				//DONE: always keep a RC_TData pointer around, with additional info about guides
				RC_TData* tdata=new RC_TData(*m, ++c_tid);
				m->uptr=tdata;
				guides_RC_tdata.Add(tdata);
				if (ballgown) { //already gather exon & intron info for all ref transcripts
					tdata->rc_addFeatures(c_exon_id, uexons, guides_RC_exons,
							c_intron_id, uintrons, guides_RC_introns);
				}
				GRefData& grefdata = refguides[m->gseq_id];
				grefdata.add(&gffr, m); //transcripts already sorted by location
			}
			if (verbose) {
				printTime(stderr);
				GMessage(" %d reference transcripts loaded.\n", gffr.gflst.Count());
			}
		}

		gseqNames=GffObj::names; //might have been populated already by gff data
		gffnames_ref(gseqNames);  //initialize the names collection if not guided
		bool havePtFeatures=false;

		/*******************************************
		 ** loading point-feature data
		 *******************************************/
		if (!ptff.is_empty()) {
			FILE* f=fopen(ptff.chars(),"r");
			if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
				ptff.chars());
			//                transcripts_only    sort by location?
			int numptf=loadPtFeatures(f, refpts); //adds to gseqNames->gseqs accordingly, populates refpts
			havePtFeatures=(numptf>0);
			fclose(f);
		}

		/*******************************************
		 ** input processing
		 *******************************************/
		GHash<int> hashread;      //read_name:pos:hit_index => readlist index
		GList<GffObj>* guides=NULL; //list of transcripts on a specific reference
		GList<GPtFeature>* refptfs=NULL; //list of point-features on a specific reference
		int currentstart=0, currentend=0;
		int ng_start=0;
		int ng_end=-1;
		int ptf_idx=0; //point-feature current index in the current (*refptfs)[]
		int ng=0;
		GStr lastref;
		bool no_ref_used=true;
		int lastref_id=-1; //last seen gseq_id
		// int ncluster=0; used it for debug purposes only

	#ifdef GFF_DEBUG
		for (int r=0;r<refguides.Count();++r) {
			GRefData& grefdata = refguides[r];
			for (int k=0;k<grefdata.rnas.Count();++k) {
				GMessage("#transcript #%d : %s (%d exons)\n", k, grefdata.rnas[k]->getID(), grefdata.rnas[k]->exons.Count());
				grefdata.rnas[k]->printGff(stderr);
			}
		}
		GMessage("GFF Debug mode, exiting...\n");
		exit(0);
	#endif

		/*******************************************
		 ** Initiating bundle parameters.
		 *******************************************/
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
	#endif

		GSamRecord* brec=NULL;	
		bool more_alns=true;
		TAlnInfo* tinfo=NULL; // for --merge
		int prev_pos=0;
		bool skipGseq=false;
		
		/*******************************************
		 ** Processing BAM files one by one.
		 *******************************************/
		for (int file_idx = 0; file_idx < bamcount; file_idx++) {
			/*
			{ //DEBUG ONLY
				fprintf(stderr, "bamreader.files.Get(file_idx): %s\n", bamreader.files.Get(file_idx).chars());
			}
			*/
			unispg_gp->ProcessSample(bamreader.files.Get(file_idx));
			GStr direction("");

			// Ouput unstranded graphs into bed file.
			if (graph_bed) {
				nodeunispgfname_unstrand = outfname_prefix + "_node_unstranded_unispg.bed";
				node_unispg_unstrand_bed = fopen(nodeunispgfname_unstrand.chars(), "w");		
				fprintf(node_unispg_unstrand_bed, "track name=Sample_"+GStr(file_idx)+"_node_unstranded_cov color=0,0,0 altColor=0,0,0\n");	
			}	

			// Ouput stranded graphs into DOT / BED file.
			for (int s=0; s<2; s++) {
				if (s == 0) {
					direction = "neg";
				} else if (s == 1) {
					direction = "pos";
				}
				dotfname = outfname_prefix + "_node_"+direction.chars()+"_"+GStr(file_idx)+"_unispg.dot";
				dot = fopen(dotfname.chars(), "w");

				dotfname_vec[s]->Add(dotfname);
				dot_vec[s]->Add(dot);

				if (graph_bed) {
					// Initializing BED files.
					nodecovfname = outfname_prefix + "_node_"+direction.chars()+"_lclg_"+GStr(file_idx)+".bed";
					edgecovfname = outfname_prefix + "_edge_"+direction.chars()+"_lclg_"+GStr(file_idx)+".bed";				
					node_cov_bed = fopen(nodecovfname.chars(), "w");				
					edge_cov_bed = fopen(edgecovfname.chars(), "w");
					fprintf(node_cov_bed, "track name=Sample_"+GStr(file_idx)+"_node_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(stderr, "track name=Sample_"+GStr(file_idx)+"_node_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(edge_cov_bed, "track name=junctions Sample_"+GStr(file_idx)+"_edge_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(stderr, "track name=junctions Sample_"+GStr(file_idx)+"_edge_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					nodelclgfname_vec[s]->Add(nodecovfname);
					edgelclgfname_vec[s]->Add(edgecovfname);
					node_lclg_bed_vec[s]->Add(node_cov_bed);
					edge_lclg_bed_vec[s]->Add(edge_cov_bed);


					nodecovfname = outfname_prefix + "_node_"+direction.chars()+"_nonovp_"+GStr(file_idx)+".bed";
					edgecovfname = outfname_prefix + "_edge_"+direction.chars()+"_nonovp_"+GStr(file_idx)+".bed";				
					node_cov_bed = fopen(nodecovfname.chars(), "w");				
					edge_cov_bed = fopen(edgecovfname.chars(), "w");
					fprintf(node_cov_bed, "track name=Sample_"+GStr(file_idx)+"_node_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(stderr, "track name=Sample_"+GStr(file_idx)+"_node_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(edge_cov_bed, "track name=junctions Sample_"+GStr(file_idx)+"_edge_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(stderr, "track name=junctions Sample_"+GStr(file_idx)+"_edge_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					nodenovpfname_vec[s]->Add(nodecovfname);
					edgenovpfname_vec[s]->Add(edgecovfname);
					node_novp_bed_vec[s]->Add(node_cov_bed);
					edge_novp_bed_vec[s]->Add(edge_cov_bed);


					nodecovfname = outfname_prefix + "_node_"+direction.chars()+"_unispg_"+GStr(file_idx)+".bed";
					edgecovfname = outfname_prefix + "_edge_"+direction.chars()+"_unispg_"+GStr(file_idx)+".bed";				
					node_cov_bed = fopen(nodecovfname.chars(), "w");				
					edge_cov_bed = fopen(edgecovfname.chars(), "w");
					fprintf(node_cov_bed, "track name=Sample_"+GStr(file_idx)+"_node_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(stderr, "track name=Sample_"+GStr(file_idx)+"_node_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(edge_cov_bed, "track name=junctions Sample_"+GStr(file_idx)+"_edge_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					fprintf(stderr, "track name=junctions Sample_"+GStr(file_idx)+"_edge_"+direction.chars()+"_cov color=255,0,0 altColor=0,0,255\n");
					nodeunispgfname_vec[s]->Add(nodecovfname);
					edgeunispgfname_vec[s]->Add(edgecovfname);
					node_unispg_bed_vec[s]->Add(node_cov_bed);
					edge_unispg_bed_vec[s]->Add(edge_cov_bed);
				}
			}

			bundle->Clear();
			brec=NULL;	
			more_alns=true;
			tinfo=NULL; // for --merge
			prev_pos=0;
			skipGseq=false;

			currentstart=0;
			currentend=0;
			ng_start=0;
			ng_end=-1;
			ptf_idx=0; //point-feature current index in the current (*refptfs)[]
			ng=0;
			bamreader.start_fidx(file_idx);

			while (more_alns) { 
				bool chr_changed=false;
				int pos=0;
				const char* refseqName=NULL;
				char xstrand=0;
				int nh=1;
				int hi=0;
				int gseq_id=lastref_id;  //current chr id
				bool new_bundle=false;
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
					refseqName=brec->refName();
					xstrand=brec->spliceStrand(); // tagged strand gets priority

					fprintf(stderr, "** Before setting strand: %c\n", xstrand);
					if(xstrand=='.' && (fr_strand || rf_strand)) { // set strand if stranded library
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

					/*
					if (xstrand=='.' && brec->exons.Count()>1) {
						no_xs++;
						continue; //skip spliced alignments lacking XS tag (e.g. HISAT alignments)
					}
					// I might still infer strand later */

					if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
					pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
					chr_changed=(lastref.is_empty() || lastref!=refseqName);
					if (chr_changed) {
						skipGseq=excludeGseqs.hasKey(refseqName);
						gseq_id=gseqNames->gseqs.addName(refseqName);
						if (guided) {
							if (gseq_id>=refseqCount) {
								if (verbose)
									GMessage("WARNING: no reference transcripts found for genomic sequence \"%s\"! (mismatched reference names?)\n",
										refseqName);
							}
							else no_ref_used=false;
						}

						if (alncounts.Count()<=gseq_id) {
							alncounts.Resize(gseq_id+1);
						}
						else if (alncounts[gseq_id]>0)
								GError("%s\nAlignments (%d) already found for %s !\n",
									ERR_BAM_SORT, alncounts[gseq_id], refseqName);
						prev_pos=0;
					}
					if (pos<prev_pos) GError("%s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
						ERR_BAM_SORT, brec->name(), brec->start,  pos, refseqName, prev_pos);
					prev_pos=pos;
					if (skipGseq) continue;
					alncounts[gseq_id]++;
					nh=brec->tag_int("NH");
					if (nh==0) nh=1;
					hi=brec->tag_int("HI");
					if (mergeMode) {
						//tinfo=new TAlnInfo(brec->name(), brec->tag_int("ZF"));
							tinfo=new TAlnInfo(brec->name(), brec->uval);
						GStr score(brec->tag_str("ZS"));
						if (!score.is_empty()) {
							GStr srest=score.split('|');
							if (!score.is_empty())
								tinfo->cov=score.asDouble();
							score=srest.split('|');
							if (!srest.is_empty())
								tinfo->fpkm=srest.asDouble();
							srest=score.split('|');
							if (!score.is_empty())
								tinfo->tpm=score.asDouble();
						}
					}

					if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
						new_bundle=true;
					}
				}
				else { //no more alignments
					more_alns=false;
					new_bundle=true; //fake a new start (end of last bundle)
				}

				/*****************************
				 * Condition to start processing reads in the previous bundle
				 *****************************/
				if (new_bundle || chr_changed) {
					hashread.Clear();
					if (bundle->readlist.Count()>0) { // process reads in previous bundle
						// (readthr, junctionthr, mintranscriptlen are globals)
						if (refptfs) { //point-features defined for this reference
							while (ptf_idx<refptfs->Count() && (int)(refptfs->Get(ptf_idx)->coord)<currentstart)
								ptf_idx++;
							//TODO: what if a PtFeature is nearby, just outside the bundle?
							while (ptf_idx<refptfs->Count() && (int)(refptfs->Get(ptf_idx)->coord)<=currentend) {
								bundle->ptfs.Add(refptfs->Get(ptf_idx)); //keep this PtFeature
								ptf_idx++;
							}
						}
						bundle->getReady(currentstart, currentend);

						if (gfasta!=NULL) { //genomic sequence data requested
							GFaSeqGet* faseq=gfasta->fetch(bundle->refseq.chars());
							if (faseq==NULL) {
								GError("Error: could not retrieve sequence data for %s!\n", bundle->refseq.chars());
							}
							bundle->gseq=faseq->copyRange(bundle->start, bundle->end, false, true);
						}
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
						processBundle_CREATE_UNISPG(bundle, unispg_gp, file_idx);
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

					if (chr_changed) {
						if (guided) {
							ng=0;
							guides=NULL;
							ng_start=0;
							ng_end=-1;
							if (refguides.Count()>gseq_id && refguides[gseq_id].rnas.Count()>0) {
								guides=&(refguides[gseq_id].rnas);
								ng=guides->Count();
							}
						}
						if (havePtFeatures) {
							ptf_idx=-1;
							//setup refptf
							refptfs=NULL;
							GRefPtData rd(gseq_id);
							int ridx=refpts.IndexOf(rd);
							if (ridx>=0) {
							refptfs=&(refpts[ridx].pfs);
							ptf_idx=0;
							}
						}
						lastref=refseqName;
						lastref_id=gseq_id;
						currentend=0;
					}

					if (!more_alns) {
						if (verbose) {
			#ifndef NOTHREADS
							GLockGuard<GFastMutex> lock(logMutex);
			#endif
							if (Num_Fragments) {
								printTime(stderr);
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
					currentstart=pos;
					currentend=brec->end;
					if (guides) { //guided and guides!=NULL
						ng_start=ng_end+1;
						while (ng_start<ng && (int)(*guides)[ng_start]->end < pos) {
							// for now, skip guides which have no overlap with current read
							ng_start++;
						}
						int ng_ovl=ng_start;
						//add all guides overlapping the current read and other guides that overlap them
						while (ng_ovl<ng && (int)(*guides)[ng_ovl]->start<=currentend) { //while guide overlap
							if (currentstart>(int)(*guides)[ng_ovl]->start)
								currentstart=(*guides)[ng_ovl]->start;
							if (currentend<(int)(*guides)[ng_ovl]->end)
								currentend=(*guides)[ng_ovl]->end;
							if (ng_ovl==ng_start && ng_ovl>0) { //first time only, we have to check back all possible transitive guide overlaps
								//char* geneid=(*guides)[ng_ovlstart]->getGeneID();
								//if (geneid==NULL) geneid=(*guides)[ng_ovlstart]->getGeneName();
								//if (geneid && !bgeneids.hasKey(geneid))
								//  bgeneids.shkAdd(geneid, &ng); //whatever pointer to int
								int g_back=ng_ovl; //start from the overlapping guide, going backwards
								int g_ovl_start=ng_ovl;
								while (g_back>ng_end+1) {
									--g_back;
									//if overlap, set g_back_start=g_back and update currentstart
									if (currentstart<=(int)(*guides)[g_back]->end) {
										g_ovl_start=g_back;
										if (currentstart>(int)(*guides)[g_back]->start)
											currentstart=(int)(*guides)[g_back]->start;
									}
								} //while checking previous guides that could be pulled in this bundle
								for (int gb=g_ovl_start;gb<=ng_ovl;++gb) {
									bundle->keepGuide((*guides)[gb],
											&guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
								}
							} //needed to check previous guides for overlaps
							else
							bundle->keepGuide((*guides)[ng_ovl],
									&guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
							ng_ovl++;
						} //while guide overlap
						ng_end=ng_ovl-1; //MUST update ng_end here, even if no overlaps were found
					} //guides present on the current chromosome
					bundle->refseq=lastref;
					bundle->start=currentstart;
					bundle->end=currentend;
				} //<---- new bundle started

				/*****************************
				 * current read extends the bundle
				 * 	this might not happen if a longer guide had already been added to the bundle
				 *****************************/
				if (currentend<(int)brec->end) {
					currentend=brec->end;
					if (guides) { //add any newly overlapping guides to bundle
						bool cend_changed;
						do {
							cend_changed=false;
							while (ng_end+1<ng && (int)(*guides)[ng_end+1]->start<=currentend) {
								++ng_end;
								//more transcripts overlapping this bundle?
								if ((int)(*guides)[ng_end]->end>=currentstart) {
									//it should really overlap the bundle
									bundle->keepGuide((*guides)[ng_end],
											&guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
									if(currentend<(int)(*guides)[ng_end]->end) {
										currentend=(*guides)[ng_end]->end;
										cend_changed=true;
									}
								}
							}
						} while (cend_changed);
					}
				} //adjusted currentend and checked for overlapping reference transcripts
				GReadAlnData alndata(brec, 0, nh, hi, tinfo);
				fprintf(stderr, "brec: %d\n", brec->mapped_len);
				fprintf(stderr, "brec: %d - %d\n", brec->start, brec->end);
				fprintf(stderr, ">>> before evalReadAln xstrand: %c\n", xstrand);
				bool ovlpguide=bundle->evalReadAln(alndata, xstrand);
				fprintf(stderr, ">>> after evalReadAln xstrand: %c\n", xstrand);

				/*****************************
				 * in eonly case consider read only if it overlaps guide
				 * 	check for overlaps with ref transcripts which may set xstrand
				 *****************************/
				// eonly: for mergeMode includes estimated coverage sum in the merged transcripts
				if(!eonly || ovlpguide) {
					if (xstrand=='+') alndata.strand=1;
					else if (xstrand=='-') alndata.strand=-1;
					//GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
					//countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
					//fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);
					processRead(currentstart, currentend, *bundle, hashread, alndata);
				}
			} //for each read alignment
			
			
			delete brec;
			bamreader.stop_fidx(file_idx);

			if (file_idx != 0) {
				unispg_gp -> Clear_no2gnode_unispg();
			// 	// Copy new_no2gnode_unispg to no2gnode_unispg
        		unispg_gp -> Copy_new_no2gnode_unispg_2_no2gnode_unispg();
				unispg_gp -> Clear_new_no2gnode_unispg();
			}
		}
	



	} else if (mode == APPLY_UNISPG) {
		/*******************************************
		 ** Process arguments.
		 *******************************************/
		GArgs args(argc, argv,
		"debug;help;version;dot=;d=;fr;rf;"
		"exclude=hvno:");
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
		 ** input processing
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
		UnispgGp_APPLY* drec=NULL;
		UnispgGp_APPLY* prev_drec=NULL;
		UnispgGp_APPLY* drec_pos=NULL;
		UnispgGp_APPLY* drec_neg=NULL;
		int prev_drec_pos_start=NULL;
		int prev_drec_pos_end=NULL;
		int prev_drec_neg_start=NULL;
		int prev_drec_neg_end=NULL;

		bool next_pos_neg = NULL;

		// Initialize the dot graphs vector.
		GPVec<UnispgGp_APPLY>* graphs_vec[2];
		for(int sno=0;sno<3;sno+=2) { // skip neutral bundles -> those shouldn't have junctions
			int s=sno/2; // adjusted strand due to ignoring neutral strand
			graphs_vec[s] = new GPVec<UnispgGp_APPLY>[1024];
		}

		unispgdotfname_pos = unispgdotfname_root + "_node_pos_0_unispg.dot";
		unispgdotfname_neg = unispgdotfname_root + "_node_neg_0_unispg.dot";
	/*
	{ //DEBUG ONLY
		fprintf(stderr, "unispgdotfname_root: %s\n", unispgdotfname_root.chars());
		fprintf(stderr, "unispgdotfname_pos: %s\n", unispgdotfname_pos.chars());
		fprintf(stderr, "unispgdotfname_neg: %s\n", unispgdotfname_neg.chars());
	}
	*/
		// unispgdotfname_root
		if (fileExists(unispgdotfname_pos.chars())==0 || fileExists(unispgdotfname_neg.chars())==0) {
			GError("%sError: the dot file does not exist!\n", USAGE);
		}

		bool dot_is_open = dotreader.start(unispgdotfname_pos); //setup and open DOT input file


		bool dot_is_open_pos = dotreader_pos.start(unispgdotfname_pos); //setup and open DOT input file

		bool dot_is_open_neg = dotreader_neg.start(unispgdotfname_neg); //setup and open DOT input file
	/*
	{ //DEBUG ONLY
		fprintf(stderr, "unispgdotfname: %d\n", dot_is_open);
		fprintf(stderr, "unispgdotfname_pos: %d\n", dot_is_open_pos);
		fprintf(stderr, "unispgdotfname_neg: %d\n", dot_is_open_neg);
		fprintf(stderr, "bamreader.files.Get(file_idx): %s\n", bamreader.files.Get(0).chars());
	}
	*/

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


  		cov_file_pos.open("/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/cov_pos.txt");
  		cov_file_neg.open("/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/cov_neg.txt");
		
		cov_file_pos_norm.open("/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/cov_pos_norm.txt");
  		cov_file_neg_norm.open("/Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/results/Brain/chr22/cov_neg_norm.txt");



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

						// UnispgGp_APPLY* copy_drec_pos = new UnispgGp_APPLY(drec_pos);
						graphs_vec[1]->Add(drec_pos);
						prev_drec_pos_start = drec_pos->get_refstart();
						prev_drec_pos_end = drec_pos->get_refend();

						next_graph_pos = true;
						next_graph_neg = false;
					} else if (drec_pos->get_refstart() > drec_neg->get_refstart()) {
						fprintf(stderr, "Adding `drec_neg`!!!!\n");
						
						// UnispgGp_APPLY* copy_drec_neg = new UnispgGp_APPLY(drec_neg);
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
							// UnispgGp_APPLY* copy_drec_pos = new UnispgGp_APPLY(drec_pos);
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
							// UnispgGp_APPLY* copy_drec_neg = new UnispgGp_APPLY(drec_neg);
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
							processBundle_APPLY_UNISPG(bundle, graphs_vec);
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
					graphs_vec[s] = new GPVec<UnispgGp_APPLY>[1024];
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

  		cov_file_pos.close();
  		cov_file_neg.close();

  		cov_file_pos_norm.close();
  		cov_file_neg_norm.close();
	}
    return 0;
}