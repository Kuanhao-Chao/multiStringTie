#include "definitions.h"
#include "unispg.h"
#include "rlink.h"
#include "GThreads.h"
#include "infer_tranx.h"

#include "stdio.h"

extern bool viral;
extern bool eonly; // parameter -e ; for mergeMode includes estimated coverage sum in the merged transcripts
extern bool longreads;

// int refseqCount=0; // number of reference sequences found in the guides file
// uint runoffdist=200;
extern double Num_Fragments; //global fragment counter (aligned pairs)
extern double Frag_Len;
extern double Cov_Sum;

extern bool NoMoreBundles;


#ifndef NOTHREADS
//single producer, multiple consumers
//main thread/program is always loading the producer
extern GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
extern GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
extern GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
extern int bundleWork; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue

//GFastMutex waitMutex;
extern GMutex waitMutex; // controls threadsWaiting (idle threads counter)

extern int threadsWaiting; // idle worker threads
extern GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle

extern GConditionVar haveClear; //will notify when bundle buf space available

extern GMutex queueMutex; //controls bundleQueue and bundles access

extern GFastMutex printMutex; //for writing the output to file

extern GFastMutex logMutex; //only when verbose - to avoid mangling the log output

extern GFastMutex bamReadingMutex;

extern GFastMutex countMutex;
#endif


void processBundleUnispg(BundleData* bundle, GPVec<UnispgGp>** graphs_vec);

void noMoreBundles();