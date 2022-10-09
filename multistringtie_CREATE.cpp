#include "global_params.h"
#include "CREATE_UNISPG/processOptions_C.h"

#include "CREATE_UNISPG/processBundle_C.h"

#include "CREATE_UNISPG/unispg_C.h"

/*******************************************
 ** CREATE_UNISPG specific parameters.
 *******************************************/
UnispgGp_CREATE* unispg_gp;
bool universal_splice_graph = false;
FILE* uinigraph_out = NULL;
GStr unigraphfname;

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

// /*****************************
//  * Declaring Ballgown related data structure
//  *   table indexes for Ballgown Raw Counts data (-B/-b option)
//  *****************************/
// GPVec<RC_TData> guides_RC_tdata(true); //raw count data or other info for all guide transcripts
// GPVec<RC_Feature> guides_RC_exons(true); //raw count data for all guide exons
// GPVec<RC_Feature> guides_RC_introns(true);//raw count data for all guide introns



void multistringtie_CREATE (int argc, char*argv[]) {
    /*******************************************
     ** Process arguments.
     *******************************************/
    GArgs args(argc, argv,
    "debug;help;version;conservative;multi;graph_bed;ref=;cram-ref=;keeptmp;rseq=;ptf=;fr;rf;"
    "exclude=zihvteuLx:n:j:s:D:G:C:S:l:m:o:a:c:f:p:g:M:Bb:A:E:F:T:");
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

    outfname_prefix = outfname.copy();		
    if (outfname.endsWith(".gtf")) {
        outfname_prefix.chomp(".gtf");
    }

    if (fileExists(outfname_prefix.chars())==0) {
        //directory does not exist, create it
        if (Gmkdir(outfname_prefix.chars()) && !fileExists(outfname_prefix.chars())) {
            GError("Error: cannot create directory %s!\n", outfname_prefix.chars());
        }
    }

    GStr outfname_node;
    GStr outfname_edge;
    outfname_node = outfname_prefix+"/node"; 
    outfname_edge = outfname_prefix+"/edge"; 

    if (fileExists(outfname_node.chars())==0) {
        //directory does not exist, create it
        if (Gmkdir(outfname_node.chars()) && !fileExists(outfname_node.chars())) {
            GError("Error: cannot create directory %s!\n", outfname_node.chars());
        }
    }
    if (fileExists(outfname_edge.chars())==0) {
        //directory does not exist, create it
        if (Gmkdir(outfname_edge.chars()) && !fileExists(outfname_edge.chars())) {
            GError("Error: cannot create directory %s!\n", outfname_edge.chars());
        }
    }


    fprintf(stderr, "outfname: %s\n", outfname.chars());
    fprintf(stderr, "out_dir: %s\n", out_dir.chars());
    fprintf(stderr, "tmp_path: %s\n", tmp_path.chars());
    fprintf(stderr, "cram_ref: %s\n", cram_ref.chars());
    fprintf(stderr, "tmpfname: %s\n", tmpfname.chars());
    fprintf(stderr, "genefname: %s\n", genefname.chars());


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
            nodeunispgfname_unstrand = outfname_node + "/unstranded_unispg.bed";
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
            dotfname = outfname_prefix + "_"+direction.chars()+"_"+GStr(file_idx)+"_unispg.dot";
            dot = fopen(dotfname.chars(), "w");

            dotfname_vec[s]->Add(dotfname);
            dot_vec[s]->Add(dot);

            if (graph_bed) {
                // Initializing BED files.
                nodecovfname = outfname_node+"/"+direction.chars()+"_lclg_"+GStr(file_idx)+".bed";
                edgecovfname = outfname_edge+"/"+direction.chars()+"_lclg_"+GStr(file_idx)+".bed";				
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


                nodecovfname = outfname_node+"/"+direction.chars()+"_nonovp_"+GStr(file_idx)+".bed";
                edgecovfname = outfname_edge+"/"+direction.chars()+"_nonovp_"+GStr(file_idx)+".bed";				
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


                nodecovfname = outfname_node+"/"+direction.chars()+"_unispg_"+GStr(file_idx)+".bed";
                edgecovfname = outfname_edge+"/"+direction.chars()+"_unispg_"+GStr(file_idx)+".bed";				
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
                            // for (int gb=g_ovl_start;gb<=ng_ovl;++gb) {
                            //     bundle->keepGuide((*guides)[gb],
                            //             &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
                            // }
                        } //needed to check previous guides for overlaps
                        else
                        // bundle->keepGuide((*guides)[ng_ovl],
                        //         &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
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
                                // bundle->keepGuide((*guides)[ng_end],
                                //         &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
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

        // clear refpts data, if loaded
        if (refpts.Count()>0)
            for (int i=0;i<refpts.Count();i++) {
                refpts[i].pfs.setFreeItem(true);
            }

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
            f_out=fopen(outfname.chars(), "w");
            if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
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
                if(istr) { // this is a transcript

                    fprintf(stderr, ">> calc_tpm: %f\n", calc_tpm);
                    fprintf(stderr, ">> calc_fpkm: %f\n", calc_fpkm);
                    fprintf(stderr, ">> tcov: %f\n", tcov);

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


        gffnames_unref(gseqNames); //deallocate names collection


#ifdef GMEMTRACE
        if(verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());
#endif

    }
}
