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
    // fprintf(stderr, "&& bamcount number: %d\n", bamcount);

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

    /*******************************************
     ** This is for reading dot file.
     *******************************************/
    bool graph_empty=false;
    bool graph_pos_empty=false;
    bool graph_neg_empty=false;

    GStr curr_graph_chr = "NULL";
    GStr prev_graph_chr; 

    bool next_graph_pos=true;
    bool next_graph_neg=true;

    int prev_pos=0;
    UnispgGp_APPLY* drec=NULL;
    int prev_drec_start=NULL;
    int prev_drec_end=NULL;

    UnispgGp_APPLY* prev_drec=NULL;
    UnispgGp_APPLY* drec_pos=NULL;
    UnispgGp_APPLY* drec_neg=NULL;
    int prev_drec_pos_start=NULL;
    int prev_drec_pos_end=NULL;
    int prev_drec_neg_start=NULL;
    int prev_drec_neg_end=NULL;

    bool next_pos_neg = NULL;

    unispgdotfname = unispgdotfname_root + "_0_unispg.dot";
    unispgdotfname_pos = unispgdotfname_root + "_pos_0_unispg.dot";
    unispgdotfname_neg = unispgdotfname_root + "_neg_0_unispg.dot";
/*
{ //DEBUG ONLY
    fprintf(stderr, "unispgdotfname_root: %s\n", unispgdotfname_root.chars());
    fprintf(stderr, "unispgdotfname_pos: %s\n", unispgdotfname_pos.chars());
    fprintf(stderr, "unispgdotfname_neg: %s\n", unispgdotfname_neg.chars());
}
*/
    if (fileExists(unispgdotfname.chars())==0) {
        GError("%sError: the dot file does not exist!\n", USAGE);
    }

    bool dot_is_open = dotreader.start(unispgdotfname); //setup and open DOT input file
    dotreader.set_samples();
    
    // bool dot_is_open_pos = dotreader_pos.start(unispgdotfname_pos); //setup and open DOT input file
    // bool dot_is_open_neg = dotreader_neg.start(unispgdotfname_neg); //setup and open DOT input file
    // dotreader_pos.set_samples();
    // dotreader_neg.set_samples();

    /*******************************************
     ** Initialize the dot graphs vector.
     *******************************************/
    UnispgGp_APPLY* unispgs;
    unispgs = new UnispgGp_APPLY();
    unispgs->ProcessSamples(dotreader.samples);
    // unispgs->ProcessSamples(dotreader_pos.samples);

/*
    for (int s=0; s<dotreader_pos.samples.Count(); s++) {
        fprintf(stderr, "After processing, the sample name: %s\n", dotreader_pos.samples[s].chars());
    }
    for (int s=0; s<dotreader_neg.samples.Count(); s++) {
        fprintf(stderr, "After processing, the sample name: %s\n", dotreader_neg.samples[s].chars());
    }

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

    /*******************************************
     ** read guiding transcripts from input gff file
     *******************************************/
    gseqNames=GffObj::names; //might have been populated already by gff data
    gffnames_ref(gseqNames);  //initialize the names collection if not guided
    bool havePtFeatures=false;


    bundle->Clear();
    brec=NULL;	
    more_alns=true;
    graph_empty=false;
    // next_graph=true;

    graph_pos_empty=false;
    graph_neg_empty=false;
    next_graph_pos=true;
    next_graph_neg=true;

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


    GStr cov_pos_fname = out_dir + "cov_pos.txt";
    cov_file_pos.open(cov_pos_fname.chars());
    GStr cov_neg_fname = out_dir + "cov_neg.txt";
    cov_file_neg.open(cov_neg_fname.chars());
    
    GStr cov_posnorm_fname = out_dir + "cov_pos_norm.txt";
    cov_file_pos_norm.open(cov_posnorm_fname.chars());
    GStr cov_negnorm_fname = out_dir + "cov_neg_norm.txt";
    cov_file_neg_norm.open(cov_negnorm_fname.chars());



    /*****************************
     * Initial condition
     *****************************/
    // drec_pos = dotreader_pos.next();
    // drec_neg = dotreader_neg.next();
    // if (drec_pos.ref) {

    // }

    while (!graph_empty) {
        if ((drec=dotreader.next())!=NULL) {
            // fprintf(stderr, "***** drec_pos: %s: %d - %d\n", drec->refseq.chars(), drec->get_refstart(), drec->get_refend());
            curr_graph_chr = drec->refseq;
        } else {
            // 'curr_graph_chr' remains the same
            graph_empty = true;
        }

        
        if (graph_empty) {
            process_unispgs = true;
        } else {
            if ( unispgs->graph_num[0] == 0 && unispgs->graph_num[1] == 0 ) {
                /*****************************
                 * There are no graphs inside unispgs yet.
                 *****************************/
                process_unispgs = false;
            } else {
                /*****************************
                 * There are at least one graph inside unispgs.
                 *****************************/
                if (drec->refseq == prev_graph_chr) {
                    if (drec->refstart < prev_drec_end) {
                        process_unispgs = false;
                    } else {
                        process_unispgs = true;
                    }
                } else {
                    process_unispgs = true;
                }
            }
        }






        /*****************************
         * The graphs_vec overlapping checking is done. There are no overlaps between pos / neg.  Now, we need to process it.
         *****************************/
        if (process_unispgs) {
            uint unispgs_start = unispgs->refstart;
            uint unispgs_end = unispgs->refend;
            uint bundle_start = unispgs->refstart;
            uint bundle_end = unispgs->refend;
            /*****************************
             * Reading BAM file & check reads overlapping with the `unispgs`
             *****************************/
            new_bundle=false;
            while (more_alns && !new_bundle) {
                /*****************************
                 * Processing the previous read!
                 *****************************/
                if (brec != NULL && brec->start > 1 && brec->end > 1) {

                    nh=brec->tag_int("NH");
                    if (nh==0) nh=1;
                    hi=brec->tag_int("HI");
                    GReadAlnData alndata(brec, 0, nh, hi, NULL);
                    bool ovlpguide=bundle->evalReadAln(alndata, xstrand);

                    // fprintf(stderr, ">> brec->start: %u\n", brec->start);
                    // fprintf(stderr, ">> brec->end: %u\n", brec->end);
                    // fprintf(stderr, ">> brec->end: %u\n", brec->end);
                    // fprintf(stderr, ">> nh: %d\n", nh);
                    // fprintf(stderr, ">> hi: %d\n", hi);

                    if (xstrand=='+') alndata.strand=1;
                    else if (xstrand=='-') alndata.strand=-1;

                    /*****************************
                     * Set the boundary to the widest based on reads.
                    *****************************/
                    if (brec->start < bundle_start) {
                        // unispgs->refstart = brec->start;
                        bundle_start = brec->start;
                    }
                    if (brec->end > bundle_end) {
                        // unispgs->refend = brec->end;
                        bundle_end = brec->end;
                    }
                    processRead((int)bundle_start, (int)bundle_end, *bundle, hashread, alndata);
                    if (new_bundle) {
                        break;	
                    }
                }
            
                total_read += 1;
                bool chr_changed=false;
                const char* refseqName=NULL;
                char xstrand=0;
                int nh=1;
                int hi=0;
                int gseq_id=lastref_id;  //current chr id
                //delete brec;
                /*****************************
                 * Read another read!
                 *****************************/
                if ((brec=bamreader.next())!=NULL) {
                    // fprintf(stderr, ">> brec->start: %u\n", brec->start);
                    // fprintf(stderr, ">> brec->end: %u\n", brec->end);
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

                    // fprintf(stderr, "\t>>>>> %s brec: %d - %d\n", refseqName, brec->start, brec->end);

                    // fprintf(stderr, "** refseqName：%s\n", refseqName);
                    processed_read += 1;
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

                    if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");

                    // if (refseqName != curr_graph_chr) {
                    //     exit(-1);
                    // }

                    // fprintf(stderr, "\t>>>>> refseqName: %s; curr_graph_chr: %s\n", refseqName, curr_graph_chr.chars());

                    chr_changed=(lastref.is_empty() || lastref!=refseqName);
                    if (chr_changed) {
                        new_bundle = true;
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
                        // fprintf(stderr, "Add new unstrand read!!\n");
                    }
                    
                    if (refseqName == curr_graph_chr) {
                        /*****************************
                         ** Step 1-2: Check whether reads are in the range of the graph.
                        *****************************/
                        // fprintf(stderr, "Process read (pre: %d - %d ;  now: %d - %d)!!!\n", pre_refstart, pre_refend, brec->start, brec->end);
                        if (brec->end <  unispgs_start) {
                            // ----------   |(s).................(e)|
                            // The read is outside the current bundle => skipped!
                            // fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
                            // continue;
                            // new_bundle = true;
                        }
                        if (brec->start < unispgs_start && brec->end == unispgs_start) {
                            // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
                            // fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
                            // continue;
                        } else if (brec->start < unispgs_start && brec->end > unispgs_start) {
                            // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
                            // fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
                        } else if (brec->start < unispgs_start && brec->end == unispgs_start) {
                            // -----|(s)---------(e)|
                            // fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
                        } else if (brec->start == unispgs_start && brec->start < unispgs_end && brec->end < unispgs_end) {
                            // |(s)----------.................(e)|   or   |(s)....----------........(e)|
                            // fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
                        } else if (brec->start == unispgs_start && brec->start < unispgs_end && brec->end == unispgs_end) {
                            // |(s)----------.................(e)|   or   |(s)....----------........(e)|
                            // fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
                        } else if (brec->start == unispgs_start && brec->start < unispgs_end && brec->end > unispgs_end) {
                            // |(s)----------.................(e)|   or   |(s)....----------........(e)|
                            // fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
                        } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end < unispgs_end) {
                            // |(s)----------.................(e)|   or   |(s)....----------........(e)|
                            // fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
                        } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end == unispgs_end) {
                            // |(s)----------.................(e)|   or   |(s)....----------........(e)|
                            // fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
                        } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end > unispgs_end) {
                            // |(s)----------.................(e)|   or   |(s)....----------........(e)|
                            // fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
                        } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end > unispgs_end) {
                            // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
                            // The overlapping with the current processing bundle.
                            // fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
                        } else if (brec->start > unispgs_start && brec->start == unispgs_end && brec->end > unispgs_end) {
                            // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
                            // The overlapping with the current processing bundle.
                            // fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
                            // continue;
                            // int overlap_current = 0;
                            // overlap_current = drec->get_refend() - brec->start + 1;

                            // int overlap_next = 0;
                            // overlap_next = brec->end - drec->get_refend() + 1;

                            // if (overlap_current > overlap_next) {
                            // 	// The read belongs to the current processing bundle.
                            // } else {
                            // 	// The read belongs to the next bundle or not belongs to any bundles.
                            // }
                        } else if (brec->start > unispgs_end) {
                            // fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
                            new_bundle = true;
                        }
                    }
                } else { //no more alignments
                    more_alns=false;
                    new_bundle=true; //fake a new start (end of last bundle)
                }

                /*****************************
                 * Check if we need to process the new bundle
                 *****************************/
                if (new_bundle) {
                    if (lastref.is_empty()) {
                    	lastref=refseqName;
                    	lastref_id=gseq_id;
                    }
                    // fprintf(stderr, "@@@ This is a new bundle\n");
                    // if (brec != NULL) fprintf(stderr, "\t>>>>> %s brec: %d - %d\n", refseqName, brec->start, brec->end);
                    bundle->refseq=lastref;
                    hashread.Clear();
                    if (bundle->readlist.Count()>0) { // process reads in previous bundle
                        // (readthr, junctionthr, mintranscriptlen are globals)
                        bundle->getReady((int)bundle_start, (int)bundle_end);

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
                        processBundle_APPLY_UNISPG(bundle, unispgs);
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
                    lastref = refseqName;
                    // currentstart=brec->start;
                    // currentend=brec->end;
                    // bundle->start= brec->start;
                    // bundle->readlist[0]->start;
                    // bundle->end= brec->end;
                }
            }

            /*****************************
             * Clean up the 'unispgs'.
             *****************************/
            // fprintf(stderr, "Ref: %s\n", bundle->refseq.chars());
            // for(int s=0;s<2;s+=1) { // skip neutral bundles -> those shouldn't have junctions
            //     // Before cleaning up, check the total number!!!
            //     fprintf(stderr, "\t@@ cleanup unispgs[%d]: %d\n", s, unispgs->graph_num[s]);
            //     for (int i=0; i<unispgs->graph_num[s]; i++) {
            //         fprintf(stderr, "\t\t bound: %u - %u \n", unispgs->no2gnode_unispg[s][i][1]->start, unispgs->no2gnode_unispg[s][i][ unispgs->node_nums[s][i]-2 ]->end);
            //     }
            // }
            unispgs->Clear();
        }


        /*****************************
         * Adding the current graph inside the unispgs
         *****************************/
        if (drec != NULL) {
            unispgs->AddUnispg(drec->s_single_dot, drec);
            // fprintf(stderr, "prev_graph_chr: %s;   drec->refseq: %s \n", prev_graph_chr.chars(), drec->refseq.chars());
            prev_graph_chr = drec->refseq;
            prev_drec_start = drec->refstart;
            prev_drec_end = drec->refend;            
        }
    }


//     /*****************************
//      * Processing Graphs & alignments: main algorithm
//      *      Iterating graphs
//      *****************************/
//     while (!graph_pos_empty || !graph_neg_empty) {
//         // After this, the current drec_pos is the newest positive graph.
//         if (next_graph_pos) {
//             prev_graph_chr_pos = curr_graph_chr_pos;
//             if ((drec=dotreader_pos.next())!=NULL) {
//                 // fprintf(stderr, "***** drec_pos: %d - %d\n", drec_pos->get_refstart(), drec_pos->get_refend());
//                 curr_graph_chr = drec->refseq;
//                 curr_graph_chr_pos = drec->refseq;
//             } else {
//                 graph_pos_empty = true;
//                 next_graph_pos = false;
//             }
//         }
//         // After this, the current drec_pos is the newest negative graph.
//         if (next_graph_neg) {
//             prev_graph_chr_neg = curr_graph_chr_neg;
//             if ((drec=dotreader_neg.next())!=NULL) {
//                 // fprintf(stderr, "***** drec_neg: %d - %d\n", drec_neg->get_refstart(), drec_neg->get_refend());
//                 curr_graph_chr = drec->refseq;
//                 curr_graph_chr_neg = drec->refseq;
//             } else {
//                 graph_neg_empty = true;
//                 next_graph_neg = false;
//             }
//         }

//         if (curr_graph_chr_pos == curr_graph_chr && curr_graph_chr_neg == curr_graph_chr) {
//             // check if overlap.
//             if (drec->refstart <= prev_drec_end) {

//             } else {

//             }

//         } else if (curr_graph_chr_pos == curr_graph_chr && curr_graph_chr_neg != curr_graph_chr) {
//             // check positive overlap
//             if (drec->refstart <= prev_drec_end) {

//             } else {
                
//             }
//         } else if (curr_graph_chr_pos != curr_graph_chr && curr_graph_chr_neg == curr_graph_chr) {
//             // check negative overlap
//             if (drec->refstart <= prev_drec_end) {

//             } else {
                
//             }
//         } else {
//             GERROR("Error: 'curr_graph_chr' has to be either 'curr_graph_chr_pos' or 'curr_graph_chr_neg'. \n");
//         }


//         prev_drec_start = drec->refstart;
//         prev_drec_end = drec->refend;
//         prev_graph_chr = curr_graph_chr;
//         // if (next_graph_pos) {
//         //     if ((drec_pos=dotreader_pos.next())!=NULL) {
//         //         // fprintf(stderr, "***** drec_pos: %d - %d\n", drec_pos->get_refstart(), drec_pos->get_refend());
//         //         curr_graph_chr_pos = drec_pos->refseq;
//         //     } else {
//         //         graph_pos_empty = true;
//         //         next_graph_pos = false;
//         //     }
//         // }
//         // // After this, the current drec_pos is the newest negative graph.
//         // if (next_graph_neg) {
//         //     if ((drec_neg=dotreader_neg.next())!=NULL) {
//         //         // fprintf(stderr, "***** drec_neg: %d - %d\n", drec_neg->get_refstart(), drec_neg->get_refend());
//         //         curr_graph_chr_neg = drec_neg->refseq;
//         //     } else {
//         //         graph_neg_empty = true;
//         //         next_graph_neg = false;
//         //     }
//         // }






//         if (curr_graph_chr_pos == curr_graph_chr_neg) {
//             /*****************************
//              * current positive & negative graphs are on the same chromosome.
//              *****************************/
//             last_graph_chr = curr_graph_chr;
//             curr_graph_chr = curr_graph_chr_pos;
//             // normal process.
//             if (drec_pos != NULL && drec_neg != NULL) {
//                 // Comparing pos & neg, and decide which one to chain.
//             } else if (drec_pos != NULL && drec_neg == NULL) {
//                 // Chain pos
//             } else if (drec_pos == NULL && drec_neg != NULL) {
//                 // Chain neg
//             } else if (drec_pos == NULL && drec_neg == NULL) {
//                 // No chaining. Process the last group.
//             }
//         } else {
//             if (curr_graph_chr_pos == curr_graph_chr) {
//                 // Keep processing positive graph until it goes to the next chromosome.
//                 // Chain pos
//             } else if (curr_graph_chr_neg == curr_graph_chr) {
//                 // Keep processing negative graph until it goes to the next chromosome.
//                 // Chain neg
//             } else if (curr_graph_chr == NULL) {
//                 GERROR("Error: The first positive unispg & negative unispg must be on the same chromosome. \n");
//             } else {
//             /*****************************
//              * (1) curr_graph_chr != curr_graph_chr_pos 
//              * (2) curr_graph_chr != curr_graph_chr_neg
//              * (3) curr_graph_chr != NULLcurrent positive & negative graphs are on the same chromosome.
//              * 
//              * Make curr positive chromosome as the current chromosome.
//              *****************************/
//                 last_graph_chr = curr_graph_chr;
//                 curr_graph_chr = curr_graph_chr_pos;
//                 // Keep processing positive graph until it goes to the next chromosome.
//                 // Chain pos
//             }
//         }
























//         /*****************************
//          * There are no graphs inside unispgs yet.
//          *****************************/
//         if ( (prev_drec_neg_start == NULL && prev_drec_neg_end == NULL && prev_drec_pos_start == NULL && prev_drec_pos_end == NULL) || (unispgs->graph_num[0] == 0 && unispgs->graph_num[1] == 0) ) {
//             // This is for beginning condition.
//             // fprintf(stderr, ">> This is for beginning condition.\n");
//             // fprintf(stderr, ">> Both `prev_drec_neg` and `prev_drec_pos` are NULL.\n");
//             if (!graph_pos_empty && !graph_neg_empty) {
//                 // fprintf(stderr, ">> drec_pos->get_refstart(): %u \n>> drec_neg->get_refstart(): %u \n", drec_pos->get_refstart(), drec_neg->get_refstart());
//                 process_unispgs = false;
//                 if (drec_pos->get_refstart() <= drec_neg->get_refstart()) {
//                     // fprintf(stderr, "Adding `drec_pos`!!!!\n");

//                     // UnispgGp_APPLY* copy_drec_pos = new UnispgGp_APPLY(drec_pos);
//                     // graphs_vec[1]->Add(drec_pos);
//                     unispgs->AddUnispg(1, drec_pos);
//                     prev_drec_pos_start = drec_pos->get_refstart();
//                     prev_drec_pos_end = drec_pos->get_refend();

//                     next_graph_pos = true;
//                     next_graph_neg = false;
//                 } else if (drec_pos->get_refstart() > drec_neg->get_refstart()) {
//                     // fprintf(stderr, "Adding `drec_neg`!!!!\n");
                    
//                     // UnispgGp_APPLY* copy_drec_neg = new UnispgGp_APPLY(drec_neg);
//                     // graphs_vec[0]->Add(drec_neg);
//                     unispgs->AddUnispg(0, drec_neg);
//                     prev_drec_neg_start = drec_neg->get_refstart();
//                     prev_drec_neg_end = drec_neg->get_refend();

//                     next_graph_pos = false;
//                     next_graph_neg = true;
//                 }
//             } else {
//                 process_unispgs = false;
//                 break;
//             }
//         } 
//         /*****************************
//          * There are graphs inside unispgs yet.
//          *****************************/
//         else {
//             // graphs_vec is not empty in the beginning => start chaning now!
//             if (prev_drec_neg_start != NULL && prev_drec_neg_end != NULL) {
//                 // fprintf(stderr, ">> pre_neg, cur_pos\n");
//                 // fprintf(stderr, ">> prev_drec_neg_start: %d; prev_drec_neg_end: %d\n", prev_drec_neg_start, prev_drec_neg_end);
//                 if (!graph_pos_empty) {
//                     bool overlap = segs_overlap(prev_drec_neg_start, prev_drec_neg_end, drec_pos->get_refstart(), drec_pos->get_refend());
//                     // fprintf(stderr, ">> overlap: %d", overlap);

//                     if (overlap) {
//                         process_unispgs = false;
//                         next_graph_pos = true;
//                         // fprintf(stderr, ">> Adding `drec_pos`: %d - %d \n", drec_pos->get_refstart(), drec_pos->get_refend());
//                         // fprintf(stderr, "\t\t>> `prev_drec_neg`: %d - %d \n", prev_drec_neg_start, prev_drec_neg_end);
//                         // UnispgGp_APPLY* copy_drec_pos = new UnispgGp_APPLY(drec_pos);
//                         // graphs_vec[1]->Add(drec_pos);
//                         unispgs->AddUnispg(1, drec_pos);
//                         prev_drec_pos_start = drec_pos->get_refstart();
//                         prev_drec_pos_end = drec_pos->get_refend();
//                         next_graph_pos = true;
//                         next_graph_neg = false;
//                         continue;
//                     } else {
//                         process_unispgs = true;
//                         next_graph_pos = false;
//                         next_graph_neg = false;
//                     }
//                 }
//             }

//             if (prev_drec_pos_start != NULL && prev_drec_pos_end != NULL) {
//                 // fprintf(stderr, ">> pre_pos, cur_neg\n");
//                 if (!graph_neg_empty) {
//                     bool overlap = segs_overlap(prev_drec_pos_start, prev_drec_pos_end, drec_neg->get_refstart(), drec_neg->get_refend());
//                     if (overlap) {
//                         process_unispgs = false;
//                         next_graph_neg = true;
//                         // fprintf(stderr, ">> Adding `drec_neg`: %d - %d \n", drec_neg->get_refstart(), drec_neg->get_refend());
//                         // fprintf(stderr, "\t\t>> `prev_drec_pos`: %d - %d \n", prev_drec_pos_start, prev_drec_pos_end);
//                         // UnispgGp_APPLY* copy_drec_neg = new UnispgGp_APPLY(drec_neg);
//                         // graphs_vec[0]->Add(drec_neg);
//                         unispgs->AddUnispg(0, drec_neg);
//                         prev_drec_neg_start = drec_neg->get_refstart();
//                         prev_drec_neg_end = drec_neg->get_refend();
//                         next_graph_pos = false;
//                         next_graph_neg = true;
//                         continue;
//                     } else {
//                         process_unispgs = true;
//                         next_graph_pos = false;
//                         next_graph_neg = false;
//                     }
//                 }
//             }
//         }


//         /*****************************
//          * This is for processing the last bundle.
//          *****************************/        
//         if (graph_pos_empty && graph_neg_empty) {
//             process_unispgs = true;
//         }

//         /*****************************
//          * The graphs_vec overlapping checking is done. There are no overlaps between pos / neg.  Now, we need to process it.
//          *****************************/
//         if (process_unispgs) {
//             int unispgs_start = unispgs->refstart;
//             int unispgs_end = unispgs->refend;
//             /*****************************
//              * Reading BAM file & check reads overlapping with the `unispgs`
//              *****************************/
//             bool read_in_unispg = true;
//             while (more_alns && read_in_unispg) {
//                 total_read += 1;
//                 bool chr_changed=false;
//                 int pos=0;
//                 const char* refseqName=NULL;
//                 char xstrand=0;
//                 int nh=1;
//                 int hi=0;
//                 int gseq_id=lastref_id;  //current chr id
//                 bool new_bundle=false;
//                 bool process_read=true;
//                 //delete brec;
//                 if ((brec=bamreader.next())!=NULL) {
//                     if (brec->isUnmapped()) continue;
//                     if (brec->start<1 || brec->mapped_len<10) {
//                         if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
//                                 brec->name(), brec->start, brec->mapped_len);
//                         continue;
//                     }
// #ifdef DBG_ALN_DATA
//                     dbg_waln(brec);
// #endif

//                     r_start = brec->start; //start<end always!
//                     r_end = brec->end;

//                     refseqName=brec->refName();
//                     xstrand=brec->spliceStrand(); // tagged strand gets priority

//                     // fprintf(stderr, "** refseqName：%s\n", refseqName);
//                     processed_read += 1;
//                     /*****************************
//                      * set strand if stranded library
//                      *****************************/
//                     if(xstrand=='.' && (fr_strand || rf_strand)) {
//                         if(brec->isPaired()) { // read is paired
//                             if(brec->pairOrder()==1) { // first read in pair
//                                 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
//                                 else xstrand='-';
//                             }
//                             else {
//                                 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
//                                 else xstrand='+';
//                             }
//                         }
//                         else {
//                             if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
//                             else xstrand='-';
//                         }
//                     }

//                     if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
//                     pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
//                     chr_changed=(lastref.is_empty() || lastref!=refseqName);
//                     if (chr_changed) {
//                         new_bundle = true;
//                         skipGseq=excludeGseqs.hasKey(refseqName);
//                         gseq_id=gseqNames->gseqs.addName(refseqName);

//                         if (alncounts.Count()<=gseq_id) {
//                             alncounts.Resize(gseq_id+1);
//                         }
//                         else if (alncounts[gseq_id]>0)
//                                 GError("%s\nAlignments (%d) already found for %s !\n",
//                                     ERR_BAM_SORT, alncounts[gseq_id], refseqName);
//                     }

//                     /*****************************
//                      * Found the overlapping between a read & the unispg
//                      *****************************/
//                     // if (xstrand == '+' || xstrand=='.') {
//                     if (xstrand == '+') {
//                         pos_strand += 1;
//                         pos_strand_nh += brec->tag_int("NH");
//                     } else if (xstrand=='-') {
//                         neg_strand += 1;
//                         neg_strand_nh += brec->tag_int("NH");
//                     } else if (xstrand=='.') {
//                         unstrand += 1;
//                         unstrand_nh += brec->tag_int("NH");
//                         // fprintf(stderr, "Add new unstrand read!!\n");
//                     }
//                     /*****************************
//                      ** Step 1-2: Check whether reads are in the range of the graph.
//                     *****************************/
//                     // fprintf(stderr, "Process read (pre: %d - %d ;  now: %d - %d)!!!\n", pre_refstart, pre_refend, brec->start, brec->end);
//                     // fprintf(stderr, "\t>>>>> brec: %d - %d\n", brec->start, brec->end);
//                     if (brec->end <  unispgs_start) {
//                         // ----------   |(s).................(e)|
//                         // The read is outside the current bundle => skipped!
//                         // fprintf(stderr, "\t** Bundle: ----------   |(s).................(e)|\n");
//                         process_read = false;
//                         // continue;
//                         // read_in_unispg = false;
//                         // new_bundle = true;
//                     }
//                     if (brec->start < unispgs_start && brec->end == unispgs_start) {
//                         // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
//                         // fprintf(stderr, "\t** Bundle: ----------|(s).................(e)|\n");
//                         // continue;
//                     } else if (brec->start < unispgs_start && brec->end > unispgs_start) {
//                         // ----------|(s).................(e)|   or   -----|(s)-----............(e)|
//                         // fprintf(stderr, "\t** Bundle: -----|(s)-----............(e)|\n");
//                     } else if (brec->start < unispgs_start && brec->end == unispgs_start) {
//                         // -----|(s)---------(e)|
//                         // fprintf(stderr, "\t** Bundle: -----|(s)-------(e)|\n");
//                     } else if (brec->start == unispgs_start && brec->start < unispgs_end && brec->end < unispgs_end) {
//                         // |(s)----------.................(e)|   or   |(s)....----------........(e)|
//                         // fprintf(stderr, "\t** Bundle: |(s)----------.................(e)| \n");
//                     } else if (brec->start == unispgs_start && brec->start < unispgs_end && brec->end == unispgs_end) {
//                         // |(s)----------.................(e)|   or   |(s)....----------........(e)|
//                         // fprintf(stderr, "\t** Bundle: |(s)----------(e)|\n");
//                     } else if (brec->start == unispgs_start && brec->start < unispgs_end && brec->end > unispgs_end) {
//                         // |(s)----------.................(e)|   or   |(s)....----------........(e)|
//                         // fprintf(stderr, "\t** Bundle: |(s)----------(e)|----\n");
//                     } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end < unispgs_end) {
//                         // |(s)----------.................(e)|   or   |(s)....----------........(e)|
//                         // fprintf(stderr, "\t** Bundle: |(s)....----------........(e)|\n");
//                     } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end == unispgs_end) {
//                         // |(s)----------.................(e)|   or   |(s)....----------........(e)|
//                         // fprintf(stderr, "\t** Bundle: |(s)....----------(e)|\n");
//                     } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end > unispgs_end) {
//                         // |(s)----------.................(e)|   or   |(s)....----------........(e)|
//                         // fprintf(stderr, "\t** Bundle: |(s)....----------(e)|-----\n");
//                     } else if (brec->start > unispgs_start && brec->start < unispgs_end && brec->end > unispgs_end) {
//                         // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
//                         // The overlapping with the current processing bundle.
//                         // fprintf(stderr, "\t** Bundle: |(s)...............------(e)|-----\n");
//                     } else if (brec->start > unispgs_start && brec->start == unispgs_end && brec->end > unispgs_end) {
//                         // |(s)...............------(e)|-----    or   |(s).................(e)|----------   
//                         // The overlapping with the current processing bundle.
//                         // fprintf(stderr, "\t** Bundle: |(s).................(e)|----------\n");
//                         // continue;
//                         // int overlap_current = 0;
//                         // overlap_current = drec->get_refend() - brec->start + 1;

//                         // int overlap_next = 0;
//                         // overlap_next = brec->end - drec->get_refend() + 1;

//                         // if (overlap_current > overlap_next) {
//                         // 	// The read belongs to the current processing bundle.
//                         // } else {
//                         // 	// The read belongs to the next bundle or not belongs to any bundles.
//                         // 	read_in_unispg = false;
//                         // }
//                     } else if (brec->start > unispgs_end) {
//                         // fprintf(stderr, "\t** Bundle: |(s).................(e)|   ----------\n");
//                         read_in_unispg = false;
//                         new_bundle = true;
//                         process_read = false;
//                     }
//                 } else { //no more alignments
//                     more_alns=false;
//                     new_bundle=true; //fake a new start (end of last bundle)
//                 }

//                 /*****************************
//                  * Process the new bundle
//                  *****************************/
//                 if (new_bundle) {
//                     if (lastref.is_empty()) {
//                     	lastref=refseqName;
//                     	lastref_id=gseq_id;
//                     }
//                     bundle->refseq=refseqName;
//                     fprintf(stderr, "bundle->refseq: %s\n", bundle->refseq.chars());
//                     // fprintf(stderr, "This is a new bundle\n");
//                     // fprintf(stderr, "bundle->readlist.Count(): %d\n", bundle->readlist.Count());
//                     hashread.Clear();
//                     if (bundle->readlist.Count()>0) { // process reads in previous bundle
//                         // (readthr, junctionthr, mintranscriptlen are globals)
//                         bundle->getReady(unispgs_start, unispgs_end);

//                         if (gfasta!=NULL) { //genomic sequence data requested
//                             GFaSeqGet* faseq=gfasta->fetch(bundle->refseq.chars());
//                             if (faseq==NULL) {
//                                 GError("Error: could not retrieve sequence data for %s!\n", bundle->refseq.chars());
//                             }
//                             bundle->gseq=faseq->copyRange(bundle->start, bundle->end, false, true);
//                         }
//             #ifndef NOTHREADS
//                         //push this in the bundle queue where it'll be picked up by the threads
//                         DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
//                         int qCount=0;
//                         queueMutex.lock();
//                         bundleQueue.Push(bundle);
//                         bundleWork |= 0x02; //set bit 1
//                         qCount=bundleQueue.Count();
//                         queueMutex.unlock();
//                         DBGPRINT2("##> bundleQueue.Count()=%d)\n", qCount);
//                         //wait for a thread to pop this bundle from the queue
//                         waitMutex.lock();
//                         DBGPRINT("##> waiting for a thread to become available..\n");
//                         while (threadsWaiting==0) {
//                             haveThreads.wait(waitMutex);
//                         }
//                         waitMutex.unlock();
//                         haveBundles.notify_one();
//                         DBGPRINT("##> waitMutex unlocked, haveBundles notified, current thread yielding\n");
//                         current_thread::yield();
//                         queueMutex.lock();
//                         DBGPRINT("##> queueMutex locked until bundleQueue.Count()==qCount\n");
//                         while (bundleQueue.Count()==qCount) {
//                             queueMutex.unlock();
//                             DBGPRINT2("##> queueMutex unlocked as bundleQueue.Count()==%d\n", qCount);
//                             haveBundles.notify_one();
//                             current_thread::yield();
//                             queueMutex.lock();
//                             DBGPRINT("##> queueMutex locked again within while loop\n");
//                         }
//                         queueMutex.unlock();
//             #else //no threads
//                         processBundle_APPLY_UNISPG(bundle, unispgs);
//             #endif
//                         // ncluster++; used it for debug purposes only
//                     } //have alignments to process
//                     else { //no read alignments in this bundle?
//             #ifndef NOTHREADS
//                         dataMutex.lock();
//                         DBGPRINT2("##> dataMutex locked for bundle #%d clearing..\n", bundle->idx);
//             #endif
//                         bundle->Clear();
//             #ifndef NOTHREADS
//                         dataClear.Push(bundle->idx);
//                         DBGPRINT2("##> dataMutex unlocking as dataClear got pushed idx #%d\n", bundle->idx);
//                         dataMutex.unlock();
//             #endif
//                     } //nothing to do with this bundle
//                     if (!more_alns) {
//                         if (verbose) {
//             #ifndef NOTHREADS
//                             GLockGuard<GFastMutex> lock(logMutex);
//             #endif
//                             if (Num_Fragments) {
//                                 // printTime(stderr);
//                                 GMessage(" %g aligned fragments found.\n", Num_Fragments);
//                             }
//                             //GMessage(" Done reading alignments.\n");
//                         }
//                         noMoreBundles();
//                         break;
//                     }
//         #ifndef NOTHREADS
//                     int new_bidx=waitForData(bundles);
//                     if (new_bidx<0) {
//                         //should never happen!
//                         GError("Error: waitForData() returned invalid bundle index(%d)!\n",new_bidx);
//                         break;
//                     }
//                     bundle=&(bundles[new_bidx]);
//         #endif
//                     lastref = refseqName;
//                     // currentstart=brec->start;
//                     // currentend=brec->end;
//                     // bundle->start= brec->start;
//                     // bundle->readlist[0]->start;
//                     // bundle->end= brec->end;
//                 }
            
//                 /*****************************
//                  * Actually start processing reads.
//                  *****************************/
//                 if (brec != NULL) {
//                     nh=brec->tag_int("NH");
//                     if (nh==0) nh=1;
//                     hi=brec->tag_int("HI");
//                     // Old conditions to stop adding reads into a bundle
//                     // if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
//                     // 	new_bundle=true;
//                     // }
//                     // fprintf(stderr, "brec: %d\n", brec->mapped_len);
//                     // if (brec->start > drec->get_refend()) {
//                     // fprintf(stderr, "brec: %d - %d\n", brec->start, brec->end);
//                     // 	new_bundle=true;
//                     // }
//                     GReadAlnData alndata(brec, 0, nh, hi, NULL);
//                     // fprintf(stderr, ">>> before evalReadAln xstrand: %c\n", xstrand);
//                     bool ovlpguide=bundle->evalReadAln(alndata, xstrand);
//                     // fprintf(stderr, ">>> after evalReadAln xstrand: %c\n", xstrand);

//                     if (xstrand=='+') alndata.strand=1;
//                     else if (xstrand=='-') alndata.strand=-1;
//                     //GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
//                     //countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
//                     //fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);
//                     // processRead(brec->start, brec->end, *bundle, hashread, alndata);
//                     processRead(unispgs_start, unispgs_end, *bundle, hashread, alndata);
                    
//                     if (new_bundle) {
//                         break;	
//                     }
//                 }
//             }

//             /*****************************
//              * Clean up the 'unispgs'.
//              *****************************/
//             for(int s=0;s<2;s+=1) { // skip neutral bundles -> those shouldn't have junctions
//                 // Before cleaning up, check the total number!!!
//                 fprintf(stderr, "\t@@ cleanup unispgs[%d]: %d\n", s, unispgs->graph_num[0]);
//                 for (int i=0; i<unispgs->graph_num[s]; i++) {
//                     fprintf(stderr, "\t\t bound: %u - %u \n", unispgs->no2gnode_unispg[s][i][1]->start, unispgs->no2gnode_unispg[s][i][ unispgs->node_nums[s][i]-2 ]->end);
//                 }
//             }
//             unispgs->Clear();
//         }

//         last_graph_chr_pos = curr_graph_chr_pos;
//         last_graph_chr_neg = curr_graph_chr_neg;
//     }



    
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

    // // clear refpts data, if loaded
    // if (refpts.Count()>0)
    //     for (int i=0;i<refpts.Count();i++) {
    //         refpts[i].pfs.setFreeItem(true);
    //     }

    fclose(f_out);
    if (c_out && c_out!=stdout) fclose(c_out);

    // if(verbose && no_xs>0)
    //     GMessage("Number spliced alignments missing the XS tag (skipped): %d\n",no_xs);

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
                // this is a transcript
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
