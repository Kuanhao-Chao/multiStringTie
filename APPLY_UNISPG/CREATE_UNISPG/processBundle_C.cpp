#include "processBundle_C.h"

void processBundle_CREATE_UNISPG(BundleData* bundle, UnispgGp_CREATE* unispg_gp, int fidx) {
	if (verbose) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(logMutex);
#endif
		printTime(stderr);
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

	infer_transcripts_CREATE_UNISPG(bundle, unispg_gp, fidx);

	// if (ballgown && bundle->rc_data) {
	// 	rc_update_exons(*(bundle->rc_data));
	// }
	if (bundle->pred.Count()>0 || ((eonly || geneabundance) && bundle->keepguides.Count()>0)) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		if(mergeMode) GeneNo=printMergeResults(bundle, GeneNo,bundle->refseq);
		else GeneNo=printResults(bundle, GeneNo, bundle->refseq);
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
		printTime(stderr);
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
