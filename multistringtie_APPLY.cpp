#include "APPLY_UNISPG/processOptions_A.h"

#include "APPLY_UNISPG/processBundle_A.h"

#include "APPLY_UNISPG/unispg_A.h"

void multistringtie_APPLY (int argc, char*argv[]) {
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
    
    // UnispgGp_APPLY* unispgs;
    // unispgs = new UnispgGp_APPLY();

    unispgdotfname_pos = unispgdotfname_root + "_pos_0_unispg.dot";
    unispgdotfname_neg = unispgdotfname_root + "_neg_0_unispg.dot";
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


    fprintf(stderr, "outfname: %s\n", outfname.chars());
    fprintf(stderr, "out_dir: %s\n", out_dir.chars());
    fprintf(stderr, "f_basename: %s\n", f_basename.chars());
    fprintf(stderr, "tmp_path: %s\n", tmp_path.chars());
    fprintf(stderr, "cram_ref: %s\n", cram_ref.chars());
    fprintf(stderr, "tmpfname: %s\n", tmpfname.chars());
    fprintf(stderr, "genefname: %s\n", genefname.chars());
    fprintf(stderr, "traindir: %s\n", traindir.chars());


    /*******************************************
     ** read guiding transcripts from input gff file
     *******************************************/
    gseqNames=GffObj::names; //might have been populated already by gff data
    gffnames_ref(gseqNames);  //initialize the names collection if not guided
    bool havePtFeatures=false;


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
    bool skipGseq=false;
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

    int pos_strand_nh = 0;
    int neg_strand_nh = 0;
    int unstrand_nh = 0;

    bool process_unispgs = false;


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
        // if ( (prev_drec_neg_start == NULL && prev_drec_neg_end == NULL && prev_drec_pos_start == NULL && prev_drec_pos_end == NULL) || (unispgs->graph_num[0] == 0 && unispgs->graph_num[1] == 0) ) {
            // This is for beginning condition.
            fprintf(stderr, ">> This is for beginning condition.\n");
            fprintf(stderr, ">> Both `prev_drec_neg` and `prev_drec_pos` are NULL.\n");
            if (more_graph_pos && more_graph_neg) {
                fprintf(stderr, ">> drec_pos->get_refstart(): %u \n>> drec_neg->get_refstart(): %u \n", drec_pos->get_refstart(), drec_neg->get_refstart());
                process_unispgs = false;
                if (drec_pos->get_refstart() <= drec_neg->get_refstart()) {
                    fprintf(stderr, "Adding `drec_pos`!!!!\n");

                    // UnispgGp_APPLY* copy_drec_pos = new UnispgGp_APPLY(drec_pos);
                    graphs_vec[1]->Add(drec_pos);
                    // unispgs->AddUnispg(1, drec_pos);
                    prev_drec_pos_start = drec_pos->get_refstart();
                    prev_drec_pos_end = drec_pos->get_refend();

                    next_graph_pos = true;
                    next_graph_neg = false;
                } else if (drec_pos->get_refstart() > drec_neg->get_refstart()) {
                    fprintf(stderr, "Adding `drec_neg`!!!!\n");
                    
                    // UnispgGp_APPLY* copy_drec_neg = new UnispgGp_APPLY(drec_neg);
                    graphs_vec[0]->Add(drec_neg);
                    // unispgs->AddUnispg(0, drec_neg);
                    prev_drec_neg_start = drec_neg->get_refstart();
                    prev_drec_neg_end = drec_neg->get_refend();

                    next_graph_pos = false;
                    next_graph_neg = true;
                }
            } else {
                process_unispgs = false;
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
                        process_unispgs = false;
                        next_graph_pos = true;
                        fprintf(stderr, ">> Adding `drec_pos`: %d - %d \n", drec_pos->get_refstart(), drec_pos->get_refend());
                        fprintf(stderr, "\t\t>> `prev_drec_neg`: %d - %d \n", prev_drec_neg_start, prev_drec_neg_end);
                        // UnispgGp_APPLY* copy_drec_pos = new UnispgGp_APPLY(drec_pos);
                        graphs_vec[1]->Add(drec_pos);
                        // unispgs->AddUnispg(1, drec_pos);
                        prev_drec_pos_start = drec_pos->get_refstart();
                        prev_drec_pos_end = drec_pos->get_refend();
                        next_graph_pos = true;
                        next_graph_neg = false;
                        continue;
                    } else {
                        process_unispgs = true;
                        next_graph_pos = false;
                        next_graph_neg = false;
                    }
                }
                fprintf(stderr, ">> process_unispgs: %d", process_unispgs);
            }

            if (prev_drec_pos_start != NULL && prev_drec_pos_end != NULL) {
                fprintf(stderr, ">> pre_pos, cur_neg\n");
                if (more_graph_neg) {
                    bool overlap = segs_overlap(prev_drec_pos_start, prev_drec_pos_end, drec_neg->get_refstart(), drec_neg->get_refend());
                    if (overlap) {
                        process_unispgs = false;
                        next_graph_neg = true;
                        fprintf(stderr, ">> Adding `drec_neg`: %d - %d \n", drec_neg->get_refstart(), drec_neg->get_refend());
                        fprintf(stderr, "\t\t>> `prev_drec_pos`: %d - %d \n", prev_drec_pos_start, prev_drec_pos_end);
                        // UnispgGp_APPLY* copy_drec_neg = new UnispgGp_APPLY(drec_neg);
                        graphs_vec[0]->Add(drec_neg);
                        // unispgs->AddUnispg(0, drec_neg);
                        prev_drec_neg_start = drec_neg->get_refstart();
                        prev_drec_neg_end = drec_neg->get_refend();
                        next_graph_pos = false;
                        next_graph_neg = true;
                        continue;
                    } else {
                        process_unispgs = true;
                        next_graph_pos = false;
                        next_graph_neg = false;
                    }
                }
            }
        }



        /*****************************
         * The graphs_vec overlapping checking is done. There are no overlaps between pos / neg.  Now, we need to process it.
         *****************************/
        if (process_unispgs) {
            int graphs_vec_start = 0;
            int graphs_vec_end = 0;

            if (graphs_vec[0]->First() == nullptr) {
                graphs_vec_start = graphs_vec[1]->First()->refstart;
            } else if (graphs_vec[1]->First() == nullptr) {
                graphs_vec_start = graphs_vec[0]->First()->refstart;
            } else {
                graphs_vec_start = (graphs_vec[0]->First()->refstart < graphs_vec[1]->First()->refstart) ? graphs_vec[0]->First()->refstart : graphs_vec[1]->First()->refstart;
            }

            // if (unispgs->graph_num[0] == 0) {
            //     graphs_vec_start = graphs_vec[1]->First()->refstart;
            // } else if (graphs_vec[1]->First() == nullptr) {
            //     graphs_vec_start = graphs_vec[0]->First()->refstart;
            // } else {
            //     graphs_vec_start = (graphs_vec[0]->First()->refstart < graphs_vec[1]->First()->refstart) ? graphs_vec[0]->First()->refstart : graphs_vec[1]->First()->refstart;
            // }

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
                    // if (lastref) {

                    // }
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

                    if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
                    pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
                    chr_changed=(lastref.is_empty() || lastref!=refseqName);
                    if (chr_changed) {
                        skipGseq=excludeGseqs.hasKey(refseqName);
                        gseq_id=gseqNames->gseqs.addName(refseqName);

                        if (alncounts.Count()<=gseq_id) {
                            alncounts.Resize(gseq_id+1);
                        }
                        else if (alncounts[gseq_id]>0)
                                GError("%s\nAlignments (%d) already found for %s !\n",
                                    ERR_BAM_SORT, alncounts[gseq_id], refseqName);
                    }

                    /*****************************
                     * Found the overlapping between a read & the unispg
                     *****************************/
                    // if (xstrand == '+' || xstrand=='.') {
                    if (xstrand == '+') {
                        pos_strand += 1;
                        pos_strand_nh += brec->tag_int("NH");
                    } else if (xstrand=='-') {
                        neg_strand += 1;
                        neg_strand_nh += brec->tag_int("NH");
                    } else if (xstrand=='.') {
                        unstrand += 1;
                        unstrand_nh += brec->tag_int("NH");
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
                    if (lastref.is_empty()) {
                    	lastref=refseqName;
                    	lastref_id=gseq_id;
                    }
                    bundle->refseq=lastref;                    
                    fprintf(stderr, "This is a new bundle\n");
                    fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
                    hashread.Clear();
                    if (bundle->readlist.Count()>0) { // process reads in previous bundle
                        // (readthr, junctionthr, mintranscriptlen are globals)
                        bundle->getReady(graphs_vec_start, graphs_vec_end);

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
    fprintf(stderr, "boundary_counter: %d\n", boundary_counter);
    fprintf(stderr, "skip_counter: %d\n\n\n", skip_counter);

    fprintf(stderr, "positve read num (nh): %d\n", pos_strand_nh);
    fprintf(stderr, "negative read num (nh): %d\n", neg_strand_nh);
    fprintf(stderr, "unstranded read num (nh): %d\n", unstrand_nh);
    
    fprintf(stderr, "skip_counter_nh: %d\n", skip_counter_nh);

    delete brec;
    bamreader.stop(); //close all BAM files
    delete drec;
    dotreader.stop(); //close all DOT files
    // dotreader_pos.stop(); //close all DOT files
    // dotreader_neg.stop(); //close all DOT files








    if (guided && no_ref_used) {
            GMessage("WARNING: no reference transcripts were found for the genomic sequences where reads were mapped!\n"
                    "Please make sure the -G annotation file uses the same naming convention for the genome sequences.\n");
    }

    delete gfasta;

#ifndef NOTHREADS
    for (int t=0;t<num_cpus;t++)
        threads[t].join();
    if (verbose) {
    printTime(stderr);
    GMessage(" All threads finished.\n");
    }
    delete[] threads;
    delete[] bundles;
#else
    if (verbose) {
        printTime(stderr);
        GMessage(" Done.\n");
    }
#endif

#ifdef DBG_ALN_DATA
    fclose(fdbgaln);
#endif

#ifdef B_DEBUG
    fclose(dbg_out);
#endif
    // if (mergeMode && guided )
    //     writeUnbundledGuides(refguides, f_out);


    // // clear refpts data, if loaded
    // if (refpts.Count()>0)
    //     for (int i=0;i<refpts.Count();i++) {
    //         refpts[i].pfs.setFreeItem(true);
    //     }

    fclose(f_out);
    if (c_out && c_out!=stdout) fclose(c_out);

    // if(verbose && no_xs>0)
    //     GMessage("Number spliced alignments missing the XS tag (skipped): %d\n",no_xs);

    if(!mergeMode) {
        if(verbose) {
            GMessage("Total count of aligned fragments: %g\n", Num_Fragments);
            if (Num_Fragments)
            GMessage("Fragment coverage length: %g\n", Frag_Len/Num_Fragments);
        }

        f_out=stdout;
        if(outfname!="stdout") {
            fprintf(stderr, "outfname: %s\n", outfname.chars());
            GStr unispg_f_out = out_dir+"UNISPG_"+f_basename;
            f_out=fopen(unispg_f_out.chars(), "w");
            if (f_out==NULL) GError("Error creating output file %s\n", unispg_f_out.chars());
        }

        fprintf(f_out,"# ");
        args.printCmdLine(f_out);
        fprintf(f_out,"# StringTie version %s\n",VERSION);

        //fprintf(stderr,"cov_sum=%g frag_len=%g num_frag=%g\n",Cov_Sum,Frag_Len,Num_Fragments);

        FILE *g_out=NULL;
        if(geneabundance) {
            g_out=fopen(genefname.chars(),"w");
            if (g_out==NULL)
                GError("Error creating gene abundance output file %s\n", genefname.chars());
            fprintf(g_out,"Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM\n");
        }

        FILE* ftmp_in=fopen(tmpfname.chars(),"rt");
        if (ftmp_in!=NULL) {
            char* linebuf=NULL;
            int linebuflen=5000;
            GMALLOC(linebuf, linebuflen);
            int nl;
            int istr;
            int tlen;
            float tcov; //do we need to increase precision here ? (double)
            float calc_fpkm;
            float calc_tpm;
            int t_id;
            while(fgetline(linebuf,linebuflen,ftmp_in)) {
                sscanf(linebuf,"%d %d %d %d %g", &istr, &nl, &tlen, &t_id, &tcov);
                if (tcov<0) tcov=0;
                if (Frag_Len>0.001) calc_fpkm=tcov*1000000000/Frag_Len;
                    else calc_fpkm=0.0;
                if (Cov_Sum>0.00001) calc_tpm=tcov*1000000/Cov_Sum;
                    else calc_tpm=0.0;
                if(istr) { 
                    fprintf(stderr, ">> calc_tpm: %f\n", calc_tpm);
                    fprintf(stderr, ">> calc_fpkm: %f\n", calc_fpkm);
                    fprintf(stderr, ">> tcov: %f\n", tcov);

                    // this is a transcript
                    // if (ballgown && t_id>0) {
                    //     guides_RC_tdata[t_id-1]->fpkm=calc_fpkm;
                    //     guides_RC_tdata[t_id-1]->cov=tcov;
                    // }
                    for(int i=0;i<nl;i++) {
                        fgetline(linebuf,linebuflen,ftmp_in);
                        if(!i) {
                            //linebuf[strlen(line)-1]='\0';
                            fprintf(f_out,"%s",linebuf);
                            fprintf(f_out," FPKM \"%.6f\";",calc_fpkm);
                            fprintf(f_out," TPM \"%.6f\";",calc_tpm);
                            fprintf(f_out,"\n");
                        }
                        else fprintf(f_out,"%s\n",linebuf);
                    }
                }
                else { // this is a gene -> different file pointer
                    fgetline(linebuf, linebuflen, ftmp_in);
                    fprintf(g_out, "%s\t%.6f\t%.6f\n", linebuf, calc_fpkm, calc_tpm);
                }
            }
            // if (guided) {
            //     writeUnbundledGuides(refguides, f_out, g_out);
            // }
            fclose(f_out);
            fclose(ftmp_in);
            if(geneabundance) fclose(g_out);
            GFREE(linebuf);
            if (!keepTempFiles) {
                remove(tmpfname.chars());
            }
        }
        else {
            fclose(f_out);
            GError("No temporary file %s present!\n",tmpfname.chars());
        }

        //lastly, for ballgown, rewrite the tdata file with updated cov and fpkm
        // if (ballgown) {
        //     rc_writeRC(guides_RC_tdata, guides_RC_exons, guides_RC_introns,
        //             f_tdata, f_edata, f_idata, f_e2t, f_i2t);
        // }
    }

    if (!keepTempFiles) {
        tmp_path.chomp('/');
        remove(tmp_path);
    }


    // gffnames_unref(gseqNames); //deallocate names collection


#ifdef GMEMTRACE
    if(verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());
#endif


















    cov_file_pos.close();
    cov_file_neg.close();

    cov_file_pos_norm.close();
    cov_file_neg_norm.close();
}
