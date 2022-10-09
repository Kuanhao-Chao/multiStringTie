#include "processTransfrags_A.h"


void process_transfrags_APPLY_UNISPG(int s, int gno,int edgeno,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag,CTreePat *tr2no, GIntHash<int> &gpos, GList<CPrediction>& pred) {

	// /*
	{ // DEBUG ONLY
		printTime(stderr);
		printTime(stderr);
		fprintf(stderr,"\n***********************************************\n");
		fprintf(stderr,"*********** Inside process_transfrags_unispg: \n");
		fprintf(stderr,"***********************************************\n");
		fprintf(stderr,"There are %d transfrags before clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d](%f,%f) long=%d short=%d usepath=%d:",i,transfrag[i]->abundance,transfrag[i]->srabund,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	// */

	GPVec<CTransfrag> srfrag(false);
	// eliminate transfrags below threshold (they represent noise) if they don't come from source
	if(!eliminate_transfrags_under_thr(gno,gpos,transfrag,tr2no,trthr,srfrag) && srfrag.Count()) { // long transfrags but only to source/sink
		srfrag.Clear();
	}

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d] abund=%f long:%d: short:%d usepath:%d",i,transfrag[i]->abundance,transfrag[i]->longread,transfrag[i]->shortread,(int)transfrag[i]->usepath);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr,"\n");
		}
	}
	*/

	GBitVec allpat(gno+edgeno);

	bool trsort=true;

	if(longreads) { // add source/sink links but only if they need to be added to explain the traversals in the graph
		transfrag.Sort(longtrCmp); // most abundant transfrag in the graph come first, then the ones with most nodes, then the ones more complete
		trsort=false;
		int source=0;
		int sink=gno-1;
		GVec<int> hassource(gno,-1); // remembers transcript number that links given node to source
		GVec<int> hassink(gno,-1); // remembers transcript number that links given node to sink
		GBitVec keepsource(gno); // if not set then I can remove link from node to source; keeps source link if it exists otherwise otherwise
		GBitVec keepsink(gno); // if not set then I can remove link from node to sink; keeps sink link if it exists otherwise otherwise
		GVec<CLongTrf> keeptrf; // keeps all potential transfrags that will be kept from most abundant to least, unassembled
		float zero=0;
		GVec<float> addsource(gno,zero);
		GVec<float> addsink(gno,zero);
		int edgedist=CHI_WIN; // I need to be consistent (if I change here then I need to change in update_abundance too)
		int ssdist=longintronanchor;

		/*
		{ // DEBUG ONLY
			//printTime(stderr);
			fprintf(stderr,"\nThere are %d transfrags after clean up:\n",transfrag.Count());
			for(int i=0;i<transfrag.Count();i++) {
				fprintf(stderr,"transfrag[%d] abund=%f:",i,transfrag[i]->abundance);
				for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
				if(transfrag[i]->guide) fprintf(stderr," guide");
				fprintf(stderr,"\n");
			}
		}
		*/

		for(int t1=0;t1<transfrag.Count();t1++) {
			/*fprintf(stderr,"Consider t=%d with abund=%f guide=%d and nodes:",t1,transfrag[t1]->abundance,transfrag[t1]->guide);
			for(int j=0;j<transfrag[t1]->nodes.Count();j++) {
				if(j) {
					int *pos=gpos[edge(transfrag[t1]->nodes[j-1],transfrag[t1]->nodes[j],gno)];
					if(pos && transfrag[t1]->pattern[*pos]) {
						fprintf(stderr,"-");
					}
					else fprintf(stderr," ");
				}
				fprintf(stderr,"%d",transfrag[t1]->nodes[j]);
			} fprintf(stderr,"\n");*/
			if(!transfrag[t1]->nodes[0]) {
				hassource[transfrag[t1]->nodes[1]]=t1;
				//fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",transfrag[t1]->nodes[1],t1,transfrag[t1]->abundance);
			}
			else if(transfrag[t1]->nodes.Last()==gno-1) {
				hassink[transfrag[t1]->nodes[0]]=t1;
				//fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",transfrag[t1]->nodes[0],t1,transfrag[t1]->abundance);
			}
			else {
				if(eonly && !transfrag[t1]->guide) continue; // do not remember transfrags that are not guides
				if(!keepsource[transfrag[t1]->nodes[0]]) {
					if(transfrag[t1]->longstart) keepsource[transfrag[t1]->nodes[0]]=1; //fprintf(stderr,"keep source %d\n",transfrag[t1]->nodes[0]);}
					else if(no2gnode[transfrag[t1]->nodes[0]]->hardstart) keepsource[transfrag[t1]->nodes[0]]=1;
				}
				if(!keepsink[transfrag[t1]->nodes.Last()]) {
					if(transfrag[t1]->longend) keepsink[transfrag[t1]->nodes.Last()]=1;//fprintf(stderr,"keep sink %d\n",transfrag[t1]->nodes.Last());}
					else if(no2gnode[transfrag[t1]->nodes.Last()]) keepsink[transfrag[t1]->nodes.Last()]=1;
				}
				bool included=false;
				// a transfrag that starts at source and ends at sink can never be included in a kept transfrag, so I am safe to do next
				for(int t2=0; t2<keeptrf.Count();t2++) {
					int t[2]={t1,keeptrf[t2].t}; // t1 current, t2 the one I kept
					int len[4]={MAX_NODE,MAX_NODE,MAX_NODE,MAX_NODE};
					int ret=compatible_long_APPLY_UNISPG(t,len,transfrag,no2gnode,gno,gpos);
					//fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
					if(ret){
						switch(ret) {
						case 1: // t[0] includes t[1]: it extends with introns on on one or both sides -> keep unless I can eliminate a previous one
							if(!transfrag[t[1]]->guide && transfrag[t1]->longstart && transfrag[t1]->longend && // t[1] might be included in t[0] so I might eliminate if it doesn't pass threshold
									(!no2gnode[transfrag[t[1]]->nodes[0]]->hardstart || transfrag[t[0]]->nodes[0] == transfrag[t[1]]->nodes[0] ) &&
									(!no2gnode[transfrag[t[1]]->nodes.Last()]->hardend || transfrag[t[0]]->nodes.Last() == transfrag[t[1]]->nodes.Last())) {
								//if(transfrag[t[0]]->abundance>(1-ERROR_PERC/DROP)*transfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
								if(transfrag[t[0]]->abundance>DROP*transfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
									if(len[1]<ssdist && len[3]<ssdist) { // prefer t[0] instead of t[1]
										keeptrf[t2].t=t1;
										keeptrf[t2].cov+=transfrag[t1]->abundance;
										keeptrf[t2].group.Add(t1);
										included=true; // I do not want to store transcript
										fprintf(stderr,"trf %d includes %d;  cov: %f\n",t[0],t[1], keeptrf[t2].cov);
									}
								}
							}
							break;
						case 2: // t[1] includes t[0]: extends with introns past ends of t[0] (t[1] possibly includes t[0]); t[1] is more abundant than t0
							//if(transfrag[t[1]]->guide || transfrag[t[1]]->abundance>(1-ERROR_PERC/DROP)*transfrag[t[0]]->abundance) {
							if(!transfrag[t[0]]->guide &&
									(!no2gnode[transfrag[t[0]]->nodes[0]]->hardstart || transfrag[t[0]]->nodes[0] == transfrag[t[1]]->nodes[0]) &&
									(!no2gnode[transfrag[t[0]]->nodes.Last()]->hardend || transfrag[t[0]]->nodes.Last() == transfrag[t[1]]->nodes.Last())) {
								if(len[1]<ssdist && len[3]<ssdist) {
									keeptrf[t2].cov+=transfrag[t1]->abundance;
									keeptrf[t2].group.Add(t1);
									included=true;
									fprintf(stderr,"trf %d %d equivalent start/ends; %f\n",t[1],t[0], keeptrf[t2].cov);
								}
							}
							//}
							break;
						case 3: // t1 and t0 are compatible --> just look for the edges; t1 goes further apart
							if(transfrag[t[0]]->nodes[0]!=transfrag[t[1]]->nodes[0] && transfrag[t[0]]->nodes.Last()!=transfrag[t[1]]->nodes.Last() &&
									((no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && !no2gnode[transfrag[t[1]]->nodes[0]]->hardstart &&
											!no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && no2gnode[transfrag[t[1]]->nodes.Last()]->hardend) ||
											(!no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[1]]->nodes[0]]->hardstart &&
													no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && !no2gnode[transfrag[t[1]]->nodes.Last()]->hardend))) {
								// these two transcripts both have one good start and one different --> keep them both (different option would be to add them to another transfrag that is compatible and has both hardends but then it's more complicated
								break;
							}

							// I keep both if both are guides
							//if((!transfrag[t[0]]->guide || !transfrag[t[1]]->guide) && abs(len[0])<edgedist && abs(len[2])<edgedist) { // close by
							if((!transfrag[t[0]]->guide || !transfrag[t[1]]->guide) && len[0]<edgedist && len[2]<edgedist) { // close by
								if(transfrag[t[0]]->guide || (!transfrag[t[1]]->guide && no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[0]]->nodes.Last()]->hardend))
									keeptrf[t2].t=t1; // t[0] to replace t[1]
								keeptrf[t2].cov+=transfrag[t1]->abundance;
								keeptrf[t2].group.Add(t1);
								included=true;
								fprintf(stderr,"trf %d %d equivalent start/ends; %f\n",t[1],t[0], keeptrf[t2].cov);
							}
							break;
						}
						if(included) break; // break from for loop
					}
				}

				if(!included){
					if(transfrag[t1]->guide || ((transfrag[t1]->longstart || no2gnode[transfrag[t1]->nodes[0]]->hardstart)&&
							(transfrag[t1]->longend || no2gnode[transfrag[t1]->nodes.Last()]->hardend))) { // if this is not included and has correct start/end
						CLongTrf kt(t1,transfrag[t1]->abundance);
						keeptrf.Add(kt);
						keeptrf.Last().group.Add(t1);
						//fprintf(stderr,"keep transfrag %d\n",t1);
					}
					else { // incomplete transcript, possibly wrong
						transfrag[t1]->weak=1;
						//fprintf(stderr,"Incomplete transcript %d\n",t1);
					}
				}
			}
		}


		char sign='-';
		if(s) { sign='+';}
		//GBitVec guidesource(gno);
		//GBitVec guidesink(gno);
		//for(int i=0;i<keeptrf.Count();i++) {
		for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier

			//fprintf(stderr,"Build source/sink for transfrag %d\n",keeptrf[i].t);
			int n1=transfrag[keeptrf[i].t]->nodes[0];
			int n2=transfrag[keeptrf[i].t]->nodes.Last();
			//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for transfrag %d\n",keeptrf[i].t);
			/******* previous implementation here

			addsource[n1]=keeptrf[i].cov;
			addsink[n2]=keeptrf[i].cov;
			 *******/

			/*
			if(!addsource[n1] && hassource[n1]<0) {
				int startpos=no2gnode[n1]->start-refstart;
				if(startpos-CHI_THR<0 || startpos+CHI_THR>bpcov->Count()) addsource[n1]=1;
				else {
					addsource[n1]=(get_cov(1,startpos,startpos+CHI_THR-1,bpcov)-get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov)-
							get_cov(1,startpos-CHI_THR,startpos-1,bpcov)+get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov))/(DROP*CHI_THR);
				}
			}

			if(!addsink[n2]  && hassink[n1]<0) {
				int endpos=no2gnode[n2]->end-refstart;
				if(endpos-CHI_THR<0 || endpos+CHI_THR>bpcov->Count()-1) addsink[n2]=1;
				else {
					addsink[n2]=(get_cov(1,endpos-CHI_THR+1,endpos,bpcov)-get_cov(2-2*s,endpos-CHI_THR+1,endpos,bpcov)-
							get_cov(1,endpos+1,endpos+CHI_THR,bpcov)+get_cov(2-2*s,endpos+1,endpos+CHI_THR,bpcov))/(DROP*CHI_THR);
				}
			}*/

			// all
			for(int j=0;j<keeptrf[i].group.Count();j++) {
				if(n1==transfrag[keeptrf[i].group[j]]->nodes[0]) {
					addsource[n1]+=transfrag[keeptrf[i].group[j]]->abundance;
					//fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],transfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
				}
				if(n2==transfrag[keeptrf[i].group[j]]->nodes.Last()) {
					addsink[n2]+=transfrag[keeptrf[i].group[j]]->abundance;
					//fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],transfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
				}
			}


			if(rawreads) {
				GVec<GSeg> exons;
				int j=0;
				int len=0;
				int t=keeptrf[i].t;
				while(j<transfrag[t]->nodes.Count()) {
					int nodestart=no2gnode[transfrag[t]->nodes[j]]->start;
					int nodeend=no2gnode[transfrag[t]->nodes[j]]->end;
					len+=nodeend-nodestart+1;
					while(j+1<transfrag[t]->nodes.Count() && no2gnode[transfrag[t]->nodes[j]]->end+1==no2gnode[transfrag[t]->nodes[j+1]]->start) {
						j++;
						len+=no2gnode[transfrag[t]->nodes[j]]->len();
						nodeend=no2gnode[transfrag[t]->nodes[j]]->end;
					}
					GSeg exon(nodestart,nodeend);
					exons.Add(exon);
					j++;
				}
				uint tstart=exons[0].start;
				uint tend=exons.Last().end;
				if(transfrag[t]->longstart>tstart) {
					len-=transfrag[t]->longstart-tstart;
					tstart=transfrag[t]->longstart;
				}
				if(transfrag[t]->longend && transfrag[t]->longend<tend) {
					len-=tend-transfrag[t]->longend;
					tend=transfrag[t]->longend;
				}

				CPrediction *p=new CPrediction(s, NULL, tstart, tend, keeptrf[i].cov, sign, len);
				exons[0].start=tstart;
				exons.Last().end=tend;
				p->exons=exons;
				pred.Add(p);
			}
		}
		for(int i=keeptrf.Count()-1;i>=0;i--) {
			int n1=transfrag[keeptrf[i].t]->nodes[0];
			int n2=transfrag[keeptrf[i].t]->nodes.Last();
		}

		/*
		{ // DEBUG ONLY
			fprintf(stderr,"%d keeptrf:\n",keeptrf.Count());
			for(int i=0;i<keeptrf.Count();i++) {
				fprintf(stderr,"(%d) %d abund=%f keepcov=%f",i,keeptrf[i].t,transfrag[keeptrf[i].t]->abundance,keeptrf[i].cov);
				for(int j=0;j<transfrag[keeptrf[i].t]->nodes.Count();j++) {
					fprintf(stderr," %d",transfrag[keeptrf[i].t]->nodes[j]);
				}
				fprintf(stderr,"\n");

			}
			fprintf(stderr,"%d trflong:",trflong.Count());
			for(int i=0;i<trflong.Count();i++)
				fprintf(stderr," %d",trflong[i]);
			fprintf(stderr,"\n");
		}
		*/


		// add source/sink connections
		for(int i=1;i<gno-1;i++) {
			//fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f\n",i,hassource[i],addsource[i],hassink[i],addsink[i]);
			if(hassource[i]<0) {
				if(addsource[i] && no2gnode[i]->hardstart) {
					/*float abund=trthr;
					if(no2gnode[i]->hardstart) { // node i doesn't have source as parent but it should
						abund=addsource[i];
					}*/
					// else I am not that confident this is a start
					no2gnode[i]->parent.Insert(0,zero);
					no2gnode[0]->child.Add(i);
					GVec<int> nodes;
					nodes.Add(source);
					nodes.Add(i);
					//CTransfrag *t=new CTransfrag(nodes,allpat,abund);
					CTransfrag *t=new CTransfrag(nodes,allpat,addsource[i]);
					t->pattern[source]=1;
					t->pattern[i]=1;
					t->longread=true;
					transfrag.Add(t);
					//fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
				}
			}
			else {
				if(!keepsource[i]) { // it was a mistake to link this to source --> remove
					/*fprintf(stderr,"delete abundance of trf:");
					for(int j=0;j<transfrag[hassource[i]]->nodes.Count();j++) fprintf(stderr," %d",transfrag[hassource[i]]->nodes[j]);fprintf(stderr,"\n");*/
					transfrag[hassource[i]]->abundance=0;
				}
				else { // this one has source and it should be kept => update trf abundance to a more realistic value
					//fprintf(stderr,"source %d:%d abund=%f addsource=%f\n",i,no2gnode[i]->start,transfrag[hassource[i]]->abundance,addsource[i]);
					//if(guidesource[i] && transfrag[hassource[i]]->abundance<addsource[i])
					transfrag[hassource[i]]->abundance=addsource[i];
				}
			}

			if(hassink[i]<0) {
				if(addsink[i] && no2gnode[i]->hardend) {
					/*float abund=trthr;
					if(no2gnode[i]->hardend) {
						abund=addsink[i];
					}*/
					// else  not very confident about this end

					no2gnode[i]->child.Add(sink);
					no2gnode[sink]->parent.Add(i);
					GVec<int> nodes;
					nodes.Add(i);
					nodes.Add(sink);
					CTransfrag *t=new CTransfrag(nodes,allpat,addsink[i]);
					t->pattern[sink]=1;
					t->pattern[i]=1;
					t->longread=true;
					transfrag.Add(t);
					//fprintf(stderr,"Add link to sink: %d-%d with abundance=%f\n",i,sink,addsink[i]);
				}
			}
			else {
				if(!keepsink[i]) { // it was a mistake to link this to sink --> remove
					/*fprintf(stderr,"delete abundance of trf:");
					for(int j=0;j<transfrag[hassink[i]]->nodes.Count();j++) fprintf(stderr," %d",transfrag[hassink[i]]->nodes[j]);fprintf(stderr,"\n");*/
					transfrag[hassink[i]]->abundance=0;
				}
				else {
					//fprintf(stderr,"sink %d:%d abund=%f addsink=%f\n",i,no2gnode[i]->end,transfrag[hassink[i]]->abundance,addsink[i]);
					//if(guidesink[i] && transfrag[hassink[i]]->abundance<addsink[i])
					transfrag[hassink[i]]->abundance=addsink[i];
				}
			}
		}
	}
	else if(srfrag.Count() && mixedMode) { // add source/sink links but only if they need to be added to explain the traversals in the graph
			srfrag.Sort(longtrCmp); // most abundant transfrags in the graph come first, then the ones with most nodes, then the ones more complete
			int source=0;
			int sink=gno-1;
			GVec<int> hassource(gno,-1); // remembers transcript number that links given node to source
			GVec<int> hassink(gno,-1); // remembers transcript number that links given node to sink
			GBitVec keepsource(gno); // if not set then I can remove link from node to source; keeps source link if it exists otherwise otherwise
			GBitVec keepsink(gno); // if not set then I can remove link from node to sink; keeps sink link if it exists otherwise otherwise
			GVec<CLongTrf> keeptrf; // keeps all potential transfrags that will be kept from most abundant to least, unassembled
			float zero=0;
			GVec<float> addsource(gno,zero);
			GVec<float> addsink(gno,zero);
			int edgedist=CHI_WIN; // I need to be consistent (if I change here then I need to change in update_abundance too)
			int ssdist=longintronanchor;
			int ntrf=0; // number of trflong transcripts

			/*
			{ // DEBUG ONLY
				//printTime(stderr);
				fprintf(stderr,"\nThere are %d longsrtransfrags after clean up:\n",srfrag.Count());
				for(int i=0;i<srfrag.Count();i++) {
					fprintf(stderr,"srfrag[%d] longstart=%d longend=%d abund=%f t=%d:",i,srfrag[i]->longstart,srfrag[i]->longend,srfrag[i]->abundance,int(srfrag[i]->usepath));
					for(int j=0;j<srfrag[i]->nodes.Count();j++) fprintf(stderr," %d",srfrag[i]->nodes[j]);
					fprintf(stderr,"\n");
				}
			}
			*/

			for(int t1=0;t1<srfrag.Count();t1++) {
				/*fprintf(stderr,"Consider t=%d with abund=%f and nodes:",t1,srfrag[t1]->abundance);
				for(int j=0;j<srfrag[t1]->nodes.Count();j++) {
					if(j) {
						int *pos=gpos[edge(srfrag[t1]->nodes[j-1],srfrag[t1]->nodes[j],gno)];
						if(pos && srfrag[t1]->pattern[*pos]) {
							fprintf(stderr,"-");
						}
						else fprintf(stderr," ");
					}
					fprintf(stderr,"%d",srfrag[t1]->nodes[j]);
				} fprintf(stderr,"\n");*/
				if(!srfrag[t1]->nodes[0]) {
					hassource[srfrag[t1]->nodes[1]]=t1;
					//fprintf(stderr,"Node %d in t=%d with cov=%f has source\n",srfrag[t1]->nodes[1],t1,srfrag[t1]->abundance);
				}
				else if(srfrag[t1]->nodes.Last()==gno-1) {
					hassink[srfrag[t1]->nodes[0]]=t1;
					//fprintf(stderr,"Node %d in t=%d with cov=%f has sink\n",srfrag[t1]->nodes[0],t1,srfrag[t1]->abundance);
				}
				else {
					if(eonly && !srfrag[t1]->guide) continue; // do not remember transfrags that are not guides
					if(!keepsource[srfrag[t1]->nodes[0]]) {
						if(srfrag[t1]->longstart) keepsource[srfrag[t1]->nodes[0]]=1;
						else if(no2gnode[srfrag[t1]->nodes[0]]->hardstart) keepsource[srfrag[t1]->nodes[0]]=1;
					}
					if(!keepsink[srfrag[t1]->nodes.Last()]) {
						if(srfrag[t1]->longend) keepsink[srfrag[t1]->nodes.Last()]=1;
						else if(no2gnode[srfrag[t1]->nodes.Last()]) keepsink[srfrag[t1]->nodes.Last()]=1;
					}
					bool included=false;
					// a transfrag that starts at source and ends at sink can never be included in a kept transfrag, so I am safe to do next
					for(int t2=0; t2<keeptrf.Count();t2++) {
						int t[2]={t1,keeptrf[t2].t}; // t1 current, t2 the one I kept
						int len[4]={MAX_NODE,MAX_NODE,MAX_NODE,MAX_NODE};
						int ret=compatible_long_APPLY_UNISPG(t,len,srfrag,no2gnode,gno,gpos);
						//fprintf(stderr,"  ret=%d t[0]=%d t[1]=%d len[0]=%d len[1]=%d len[2]=%d len[3]=%d\n",ret,t[0],t[1],len[0],len[1],len[2],len[3]);
						if(ret){
							switch(ret) {
							case 1: // t[0] includes t[1]: it extends with introns on on one or both sides -> keep unless I can eliminate a previous one
								if(!srfrag[t[1]]->guide && srfrag[t1]->longstart && srfrag[t1]->longend && // t[1] might be included in t[0] so I might eliminate if it doesn't pass threshold
										(!no2gnode[srfrag[t[1]]->nodes[0]]->hardstart || srfrag[t[0]]->nodes[0] == srfrag[t[1]]->nodes[0] ) &&
										(!no2gnode[srfrag[t[1]]->nodes.Last()]->hardend || srfrag[t[0]]->nodes.Last() == srfrag[t[1]]->nodes.Last())) {
									//if(srfrag[t[0]]->abundance>(1-ERROR_PERC/DROP)*srfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
									if(srfrag[t[0]]->abundance>DROP*srfrag[t[1]]->abundance) { // t[0] is within limits of t[1]
										if(len[1]<ssdist && len[3]<ssdist) { // prefer t[0] instead of t[1]
											keeptrf[t2].t=t1;
											keeptrf[t2].cov+=srfrag[t1]->abundance;
											keeptrf[t2].group.Add(t1);
											included=true; // I do not want to store transcript
											//fprintf(stderr,"trf %d includes %d\n",t[0],t[1]);
										}
									}
								}
								break;
							case 2: // t[1] includes t[0]: extends with introns past ends of t[0] (t[1] possibly includes t[0]); t[1] is more abundant than t0
								//if(srfrag[t[1]]->guide || srfrag[t[1]]->abundance>(1-ERROR_PERC/DROP)*srfrag[t[0]]->abundance) {
								if(!srfrag[t[0]]->guide &&
										(!no2gnode[srfrag[t[0]]->nodes[0]]->hardstart || srfrag[t[0]]->nodes[0] == srfrag[t[1]]->nodes[0]) &&
										(!no2gnode[srfrag[t[0]]->nodes.Last()]->hardend || srfrag[t[0]]->nodes.Last() == srfrag[t[1]]->nodes.Last())) {
									if(len[1]<ssdist && len[3]<ssdist) {
										keeptrf[t2].cov+=srfrag[t1]->abundance;
										keeptrf[t2].group.Add(t1);
										included=true;
										//fprintf(stderr,"trf %d is intronic including %d\n",t[1],t[0]);
									}
								}
								//}
								break;
							case 3: // t1 and t0 are compatible --> just look for the edges; t1 goes further apart
								if(srfrag[t[0]]->nodes[0]!=srfrag[t[1]]->nodes[0] && srfrag[t[0]]->nodes.Last()!=srfrag[t[1]]->nodes.Last() &&
										((no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && !no2gnode[srfrag[t[1]]->nodes[0]]->hardstart &&
												!no2gnode[srfrag[t[0]]->nodes.Last()]->hardend && no2gnode[srfrag[t[1]]->nodes.Last()]->hardend) ||
												(!no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && no2gnode[srfrag[t[1]]->nodes[0]]->hardstart &&
														no2gnode[srfrag[t[0]]->nodes.Last()]->hardend && !no2gnode[srfrag[t[1]]->nodes.Last()]->hardend))) {
									// these two transcripts both have one good start and one different --> keep them both (different option would be to add them to another srfrag that is compatible and has both hardends but then it's more complicated
									break;
								}

								// I keep both if both are guides
								//if((!srfrag[t[0]]->guide || !srfrag[t[1]]->guide) && abs(len[0])<edgedist && abs(len[2])<edgedist) { // close by
								if((!srfrag[t[0]]->guide || !srfrag[t[1]]->guide) && len[0]<edgedist && len[2]<edgedist) { // close by
									if(srfrag[t[0]]->guide || (!srfrag[t[1]]->guide && no2gnode[srfrag[t[0]]->nodes[0]]->hardstart && no2gnode[srfrag[t[0]]->nodes.Last()]->hardend))
										keeptrf[t2].t=t1; // t[0] to replace t[1]
									keeptrf[t2].cov+=srfrag[t1]->abundance;
									keeptrf[t2].group.Add(t1);
									included=true;
									//fprintf(stderr,"trf %d %d equivalent start/ends\n",t[1],t[0]);
								}
								break;
							}
							if(included) break; // break from for loop
						}
					}

					if(!included){
						if(srfrag[t1]->guide || ((srfrag[t1]->longstart || no2gnode[srfrag[t1]->nodes[0]]->hardstart)&&
								(srfrag[t1]->longend || no2gnode[srfrag[t1]->nodes.Last()]->hardend))) { // if this is not included and has correct start/end
							CLongTrf kt(t1,srfrag[t1]->abundance);
							keeptrf.Add(kt);
							keeptrf.Last().group.Add(t1);
							//fprintf(stderr,"keep srfrag %d\n",t1);
						}
						else { // incomplete transcript, possibly wrong
							srfrag[t1]->weak=1;
							//fprintf(stderr,"Incomplete transcript %d\n",t1);
						}
					}
				}
			}


			//GBitVec guidesource(gno);
			//GBitVec guidesink(gno);
			//for(int i=0;i<keeptrf.Count();i++) {
			for(int i=keeptrf.Count()-1;i>=0;i--) { // I add the kept transcripts to trflong from least significant to most in order to make deletion easier

				//fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);
				int n1=srfrag[keeptrf[i].t]->nodes[0];
				int n2=srfrag[keeptrf[i].t]->nodes.Last();
				//if(hassource[n1]<0 || hassink[n2]<0) fprintf(stderr,"Build source/sink for srfrag %d\n",keeptrf[i].t);

				if(!srfrag[keeptrf[i].t]->guide && ((!no2gnode[n1]->hardstart && (hassource[n1]<0 || !keepsource[n1])) ||
						(!no2gnode[n2]->hardend && (hassink[n2]<0 || !keepsink[n2])))) {
				//if((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) && (no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2]))) {
					srfrag[keeptrf[i].t]->usepath=-2-ntrf; // order of transfrag
					//fprintf(stderr,"keeptrf[%d].t=%d srfrag.usepath=%d\n",i,keeptrf[i].t,t);
					ntrf++;
					//trflong.Add(keeptrf[i].t);
				}

				/******* previous implementation here

				addsource[n1]=keeptrf[i].cov;
				addsink[n2]=keeptrf[i].cov;
				 *******/

				/*
				if(!addsource[n1] && hassource[n1]<0) {
					int startpos=no2gnode[n1]->start-refstart;
					if(startpos-CHI_THR<0 || startpos+CHI_THR>bpcov->Count()) addsource[n1]=1;
					else {
						addsource[n1]=(get_cov(1,startpos,startpos+CHI_THR-1,bpcov)-get_cov(2-2*s,startpos,startpos+CHI_THR-1,bpcov)-
								get_cov(1,startpos-CHI_THR,startpos-1,bpcov)+get_cov(2-2*s,startpos-CHI_THR,startpos-1,bpcov))/(DROP*CHI_THR);
					}
				}

				if(!addsink[n2]  && hassink[n1]<0) {
					int endpos=no2gnode[n2]->end-refstart;
					if(endpos-CHI_THR<0 || endpos+CHI_THR>bpcov->Count()-1) addsink[n2]=1;
					else {
						addsink[n2]=(get_cov(1,endpos-CHI_THR+1,endpos,bpcov)-get_cov(2-2*s,endpos-CHI_THR+1,endpos,bpcov)-
								get_cov(1,endpos+1,endpos+CHI_THR,bpcov)+get_cov(2-2*s,endpos+1,endpos+CHI_THR,bpcov))/(DROP*CHI_THR);
					}
				}*/

				// all
				for(int j=0;j<keeptrf[i].group.Count();j++) {
					if(n1==srfrag[keeptrf[i].group[j]]->nodes[0]) {
						addsource[n1]+=srfrag[keeptrf[i].group[j]]->abundance;
						//fprintf(stderr,"Add source t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n1,addsource[n1]);
					}
					if(n2==srfrag[keeptrf[i].group[j]]->nodes.Last()) {
						addsink[n2]+=srfrag[keeptrf[i].group[j]]->abundance;
						//fprintf(stderr,"Add sink t[%d]->cov=%f to node %d = %f\n",keeptrf[i].group[j],srfrag[keeptrf[i].group[j]]->abundance,n2,addsink[n2]);
						}
					}

			}
			for(int i=keeptrf.Count()-1;i>=0;i--) {
				int n1=srfrag[keeptrf[i].t]->nodes[0];
				int n2=srfrag[keeptrf[i].t]->nodes.Last();

				if(srfrag[keeptrf[i].t]->guide || ((no2gnode[n1]->hardstart || (hassource[n1]>=0 && keepsource[n1])) &&
						(no2gnode[n2]->hardend ||(hassink[n2]>=0 && keepsink[n2])))) {
					//trflong.Add(keeptrf[i].t);
					srfrag[keeptrf[i].t]->usepath=-2-ntrf; // mark transfrag that it needs to be part of trflong
					ntrf++;
				}
			}

			/*
			{ // DEBUG ONLY
				fprintf(stderr,"%d keeptrf:\n",keeptrf.Count());
				for(int i=0;i<keeptrf.Count();i++) {
					fprintf(stderr,"(%d) %d abund=%f keepcov=%f",i,keeptrf[i].t,srfrag[keeptrf[i].t]->abundance,keeptrf[i].cov);
					for(int j=0;j<srfrag[keeptrf[i].t]->nodes.Count();j++) {
						fprintf(stderr," %d",srfrag[keeptrf[i].t]->nodes[j]);
					}
					fprintf(stderr,"\n");

				}
				fprintf(stderr,"%d trflong:",trflong.Count());
				for(int i=0;i<trflong.Count();i++)
					fprintf(stderr," %d",trflong[i]);
				fprintf(stderr,"\n");
			}
			*/


			// add source/sink connections
			for(int i=1;i<gno-1;i++) {
				//fprintf(stderr,"i=%d hassource=%d addsource=%f hassink=%d addsink=%f hardstart=%d hardend=%d\n",i,hassource[i],addsource[i],hassink[i],addsink[i],no2gnode[i]->hardstart,no2gnode[i]->hardend);
				if(hassource[i]<0) {
					if(addsource[i] && no2gnode[i]->hardstart) {
						/*float abund=trthr;
						if(no2gnode[i]->hardstart) { // node i doesn't have source as parent but it should
							abund=addsource[i];
						}*/
						// else I am not that confident this is a start
						no2gnode[i]->parent.Insert(0,zero);
						no2gnode[0]->child.Add(i);
						GVec<int> nodes;
						nodes.Add(source);
						nodes.Add(i);
						//CTransfrag *t=new CTransfrag(nodes,allpat,abund);
						CTransfrag *t=new CTransfrag(nodes,allpat,addsource[i]);
						t->pattern[source]=1;
						t->pattern[i]=1;
						t->longread=true; /// this is only true for longreads;
						transfrag.Add(t);
						//fprintf(stderr,"Add link to source: 0-%d with abundance=%f\n",i, addsource[i]);
					}
				}
				else {
					/*if(keepsource[i]) {
						srfrag[hassource[i]]->abundance+=addsource[i];
					}*/
					if(!keepsource[i]) { // it was a mistake to link this to source --> remove
						/*fprintf(stderr,"delete abundance of trf:");
						for(int j=0;j<srfrag[hassource[i]]->nodes.Count();j++) fprintf(stderr," %d",srfrag[hassource[i]]->nodes[j]);fprintf(stderr,"\n");*/
						srfrag[hassource[i]]->abundance=0;
					}
					else { // this one has source and it should be kept => update trf abundance to a more realistic value
						//fprintf(stderr,"source %d:%d abund=%f addsource=%f\n",i,no2gnode[i]->start,srfrag[hassource[i]]->abundance,addsource[i]);
						//if(guidesource[i] && srfrag[hassource[i]]->abundance<addsource[i])
						srfrag[hassource[i]]->abundance=addsource[i];
					}
				}

				if(hassink[i]<0) {
					if(addsink[i] && no2gnode[i]->hardend) {
						/*float abund=trthr;
						if(no2gnode[i]->hardend) {
							abund=addsink[i];
						}*/
						// else  not very confident about this end

						no2gnode[i]->child.Add(sink);
						no2gnode[sink]->parent.Add(i);
						GVec<int> nodes;
						nodes.Add(i);
						nodes.Add(sink);
						CTransfrag *t=new CTransfrag(nodes,allpat,addsink[i]);
						t->pattern[sink]=1;
						t->pattern[i]=1;
						t->longread=true; //// this is only true for longreads
						transfrag.Add(t);
						//fprintf(stderr,"Add link to sink: %d-%d with abundance=%f\n",i,sink,addsink[i]);
					}
				}
				else {
					if(!keepsink[i]) { // it was a mistake to link this to sink --> remove
						/*fprintf(stderr,"delete abundance of trf:");
						for(int j=0;j<srfrag[hassink[i]]->nodes.Count();j++) fprintf(stderr," %d",srfrag[hassink[i]]->nodes[j]);fprintf(stderr,"\n");*/
						srfrag[hassink[i]]->abundance=0;
					}
					else {
						//fprintf(stderr,"sink %d:%d abund=%f addsink=%f\n",i,no2gnode[i]->end,srfrag[hassink[i]]->abundance,addsink[i]);
						//if(guidesink[i] && srfrag[hassink[i]]->abundance<addsink[i])
						srfrag[hassink[i]]->abundance=addsink[i];
					}
					/*if(keepsink[i]) {
						//fprintf(stderr,"Add %f to longtransfrag=%d\n",addsink[i],hassink[i]);
						srfrag[hassink[i]]->abundance+=addsink[i];
					}*/
				}
			}
	}


	// add edges between disconnected parent-child nodes
	for(int t=0;t<transfrag.Count();t++) allpat=allpat | transfrag[t]->pattern;

	for(int i=1;i<gno-1;i++) { // for all nodes check if there is a connection to child
		CGraphnodeUnispg *n=no2gnode[i];
		for(int c=0;c<n->child.Count();c++) {
			int *pos=gpos[edge(i,n->child[c],gno)];
			if(pos && !allpat[*pos]) {
				GVec<int> nodes;
				nodes.Add(i);
				nodes.Add(n->child[c]);
				GBitVec trpat(gno+edgeno);
				trpat[i]=1;
				trpat[n->child[c]]=1;
				trpat[*pos]=1;
				CTransfrag *t=new CTransfrag(nodes,trpat,trthr);
				if(longreads) t->longread=true;
				transfrag.Add(t);
			}
		}
	}

	// sort transfrag with smallest being the one that has the most nodes, and ties are decided by the abundance (largest abundance first); last transfrags all have 1 node
	if(trsort)
		transfrag.Sort(trCmp);

	// /*
	{ // DEBUG ONLY
		printTime(stderr);
		fprintf(stderr,"There are %d transfrags that remained\n",transfrag.Count());
	}
	// */

	/*
	{ // DEBUG ONLY
		//printTime(stderr);
		fprintf(stderr,"\nThere are %d transfrags after sorting:\n",transfrag.Count());
		for(int i=0;i<transfrag.Count();i++) {
			fprintf(stderr,"transfrag[%d]:",i);
			for(int j=0;j<transfrag[i]->nodes.Count();j++) fprintf(stderr," %d",transfrag[i]->nodes[j]);
			fprintf(stderr," abundance=%f usepath=%.1f\n",transfrag[i]->abundance,transfrag[i]->usepath);
		}
	}
	*/


	GVec<int> incompletetrf; //remembers incomplete transfrags (the ones that don't have edges between two consecutive nodes


	// create compatibilities
	for(int t1=0;t1<transfrag.Count();t1++) { // transfrags are processed in increasing order -> important for the later considerations

		// update nodes
		int n1=transfrag[t1]->nodes.Count();

		// fprintf(stderr, "  t1: %d;  n1: %d\n", t1, n1);
		// printBitVec(transfrag[t1]->pattern);
		// fprintf(stderr, "\n");
		if(mixedMode) {
			//transfrag[t1]->usepath=-1; // restore order
			transfrag[t1]->usepath=transfrag[t1]->abundance; // store abundance for later
		}

		if(n1>1) { // add transfrag to nodes' in and out; if a transfrag only has one node then it is not added to a node; I might want to change this for the computation of fpkm's
			bool incomplete = false;
			bool nosplice=true; // try to give less priority to unspliced reads vs spliced reads
			for(int n=0;n<n1;n++) { // for all nodes in transfrag
				if(nosplice && n) { // reduce abundance of continuous transfrags
					if(transfrag[t1]->nodes[n]!=1+transfrag[t1]->nodes[n-1] || no2gnode[transfrag[t1]->nodes[n]]->start-1!=no2gnode[transfrag[t1]->nodes[n-1]]->end) {
						nosplice=false;
					}
				}

				// fprintf(stderr, "transfrag[t1]->nodes[n-1]: %d;  transfrag[t1]->nodes[n]: %d \n", transfrag[t1]->nodes[n-1], transfrag[t1]->nodes[n]);
				if(n && n<transfrag[t1]->nodes.Count()-1) {// not first or last node
					// add t1 to in and out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// fprintf(stderr, "(1) transfrag[t1]->nodes[n]: %d\n", transfrag[t1]->nodes[n]);


					if(transfrag[t1]->nodes[n-1] && transfrag[t1]->nodes[n]<gno-1) {
						// check if transfrag t1 is incomplete between node[n-1] and node [n]
						int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];

						// fprintf(stderr, ">>> t1: %d;  pos: %d;  \n", t1, *pos);
						if(!pos || !transfrag[t1]->pattern[*pos]) {
							// there is no edge between node[n-1] and node[n]
							incomplete = assign_incomplete_trf_to_nodes_APPLY_UNISPG(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
						}
					}
				}
				else if(n) { // last but not first node
					// add t1 to in of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);

					// fprintf(stderr, "(2) transfrag[t1]->nodes[n]: %d\n", transfrag[t1]->nodes[n]);
					// fprintf(stderr, "(2) transfrag[t1]->nodes[n-1]: %d\n", transfrag[t1]->nodes[n-1]);
					// fprintf(stderr, "(2) gno %d\n", gno);

					if(transfrag[t1]->nodes[n-1] && transfrag[t1]->nodes[n]<gno-1) {
						// check if transfrag t1 is incomplete between node[n-1] and node [n]
						int *pos=gpos[edge(transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],gno)];

						// fprintf(stderr, ">>> s: %d;  gno: %d\n", s, gno);

						// fprintf(stderr, "### @@ Adding edge (%d - %d);  g_node_num: %d\n", transfrag[t1]->nodes[n-1], transfrag[t1]->nodes[n], gno);

						// fprintf(stderr, ">>> t1: %d;  \n", t1);

						// fprintf(stderr, ">>> t1: %d;  pos: %d;  \n", t1, *pos);
						if(!pos || !transfrag[t1]->pattern[*pos]) {// there is no edge between node[n-1] and node[n]
							incomplete = assign_incomplete_trf_to_nodes_APPLY_UNISPG(t1,transfrag[t1]->nodes[n-1],transfrag[t1]->nodes[n],no2gnode) or incomplete; 	// this ensures that I still have compatibilities when going through a certain path: !!! THIS IS NEW COMPARED TO PERL
						}
					}
				}
				else { // first node -> only add transfrag to out of node
					no2gnode[transfrag[t1]->nodes[n]]->trf.Add(t1);
				}
			}

			//if(nosplice) transfrag[t1]->abundance*=(1-isofrac);
			//transfrag[t1]->abundance*=0.5;

			fprintf(stderr, ">> there is no edge between node[n-1] and node[n]\n");

			if(incomplete) incompletetrf.Add(t1);
			else transfrag[t1]->real=true;

		}
		//else if(longreads) no2gnode[n1]->trf.Add(t1);
		/*
		else { // this transcript is included completely in node
			no2gnode[transfrag[t1]->nodes[0]]->frag+=transfrag[t1]->abundance;
		}
		*/
	} // end for(int t1=0;t1<transfrag.Count();t1++)

	fprintf(stderr, "After creating compatibilities!!\n");

	if(srfrag.Count() && !mixedMode) {
		srfrag.Sort(trCmp); // always start with largest super-read to solve
		for(int u=0;u<srfrag.Count();u++) //process_srfrag(srfrag[u],transfrag,no2gnode,gno,gpos);
		  if(!srfrag[u]->abundance) srfrag[u]->abundance=srfrag[u]->srabund*ERROR_PERC;
	}


	// set source-to-child transfrag abundances: optional in order not to keep these abundances too low:
	// update the abundances of the transfrags coming in from source and going to a node that doesn't have other parents than source
	// * this part was removed to improve performance
	CGraphnodeUnispg *source=no2gnode[0];
	for(int i=0;i<source->child.Count();i++) {
		float abundance=0;
		int t0=-1;
		if(no2gnode[source->child[i]]->parent.Count()==1 && !no2gnode[source->child[i]]->parent[0]) { // source is the only parent of node
			for(int j=0;j<no2gnode[source->child[i]]->trf.Count();j++) {
				int t=no2gnode[source->child[i]]->trf[j];
				if(transfrag[t]->nodes.Last()==source->child[i]) t0=t;
				else abundance+=transfrag[t]->abundance;
			}
			if(t0>-1 && transfrag[t0]->abundance) { // found transfrag from source to node and the transfrag wasn't deleted
				transfrag[t0]->abundance=abundance;
			}
		}
	}
	// */

	for(int t=0;t<incompletetrf.Count();t++)
		transfrag[incompletetrf[t]]->real=trf_real_APPLY_UNISPG(incompletetrf[t],no2gnode,transfrag,gpos,gno);


	fprintf(stderr, "Done!!\n");
}



int compatible_long_APPLY_UNISPG(int* t,int *len,GPVec<CTransfrag>& transfrag,GPVec<CGraphnodeUnispg>& no2gnode,int gno,GIntHash<int> &gpos) {
// return 0 if no compatibility, then 4*(rets)+rete where rets/rete = 1 if t[0] has extra intron, 2 if t[1] has extra intron, 3 if compatible starts/ends
// len[0] is how much t[0] extends past t[1] if positive, otherwise how much t[1] extends past t[0] at the start of transcripts
// len[1] positive: if t[0] has extra intron at start how much is t[1] inside of the intron -> default is 0; negative: same for t[1]
// len[2] is how much t[0] extends past t[1] if positive, otherwise how much t[1] extends past t[0]	at the end of transcripts
// len[3] positive: if t[0] has extra intron at end how much is t[1] inside of the intron -> default is 0; negative: same for t[1]


	if(transfrag[t[0]]->nodes[0]!=transfrag[t[1]]->nodes[0] &&
			no2gnode[transfrag[t[0]]->nodes[0]]->hardstart && no2gnode[transfrag[t[1]]->nodes[0]]->hardstart)
		return 0; // both starts can't be hard

	if(transfrag[t[0]]->nodes.Last()!=transfrag[t[1]]->nodes.Last() &&
			no2gnode[transfrag[t[0]]->nodes.Last()]->hardend && no2gnode[transfrag[t[1]]->nodes.Last()]->hardend)
		return 0; // both ends can't be hard

	uint tstart[2]={no2gnode[transfrag[t[0]]->nodes[0]]->start,no2gnode[transfrag[t[1]]->nodes[0]]->start};
	if(transfrag[t[0]]->longstart) tstart[0]=transfrag[t[0]]->longstart;
	if(transfrag[t[1]]->longstart) tstart[1]=transfrag[t[1]]->longstart;
	uint tend[2]={no2gnode[transfrag[t[0]]->nodes.Last()]->end,no2gnode[transfrag[t[1]]->nodes.Last()]->end};
	if(transfrag[t[0]]->longend) tend[0]=transfrag[t[0]]->longend;
	if(transfrag[t[1]]->longend) tend[1]=transfrag[t[1]]->longend;

	int f=0; // first starting transcript
	int s=1; // second starting transcript
	if(tstart[1]<tstart[0]) { // transfrag t1 starts before trf t0
		f=1;
		s=0;
	}

	if(transfrag[t[f]]->nodes.Last()<transfrag[t[s]]->nodes[0]) // no overlapping transcripts
		return 0;

	// check START
	int rets=3; // compatible starts by default
	int i[2]; // starting index in transfrags where they share a node
	i[0]=0;i[1]=0;
	len[0]=0;
	len[1]=0;
	while(transfrag[t[f]]->nodes[i[f]]<transfrag[t[s]]->nodes[0]) { // i[f]<transfrag[t[f]]->nodes.Count() always true due to condition of overlap before
		if(no2gnode[transfrag[t[f]]->nodes[i[f]]]->end+1<no2gnode[transfrag[t[f]]->nodes[i[f]+1]]->start) { // found intron in f
			rets=1+f;
			if(no2gnode[transfrag[t[s]]->nodes[0]]->hardstart) return 0; // there is no compatibility here
		}

		// len[0]+=no2gnode[transfrag[t[f]]->nodes[i[f]]]->len(); // this only approximates the length
		//  adj if I want real distances:
		if(i[f]) len[0]+=no2gnode[transfrag[t[f]]->nodes[i[f]]]->len();
		else len[0]+=no2gnode[transfrag[t[f]]->nodes[i[f]]]->end-tstart[f]+1;
		i[f]++;
	}

	if(transfrag[t[s]]->nodes.Last()<transfrag[t[f]]->nodes[i[f]]) return 0; // no overlap: t[s] is contained in t[f]

	// adj for real distance
	if(!i[f]) { // transfrag[t[f]]->nodes[0]==transfrag[t[s]]->nodes[0]
		len[0]=tstart[s]-tstart[f];
	}

	if(transfrag[t[s]]->nodes[0]<transfrag[t[f]]->nodes[i[f]]) { // i[f]>0 because f has to start first
		if(!no2gnode[transfrag[t[s]]->nodes[0]]->childpat[transfrag[t[f]]->nodes[i[f]]]) return 0; // I can not reach node i[f] from t[s][0]
		// to approximate:
		// len[1]=no2gnode[transfrag[t[f]]->nodes[i[f]]]->start-no2gnode[transfrag[t[s]]->nodes[0]]->start;
		// adj if I want real distances:
		len[1]=no2gnode[transfrag[t[f]]->nodes[i[f]]]->start-tstart[s];
		while(transfrag[t[s]]->nodes[i[s]]<transfrag[t[f]]->nodes[i[f]]) {
			i[s]++;
			if(no2gnode[transfrag[t[s]]->nodes[i[s]-1]]->end+1<no2gnode[transfrag[t[s]]->nodes[i[s]]]->start) { // not continuous
				return 0;
			}
		}
	}
	else if(i[f]){ // no intron -> adj lentgh for real distance
		len[0]+=tstart[s]-no2gnode[transfrag[t[s]]->nodes[0]]->start;
	}

	if(f) len[0]=-len[0]; // a negative value signals that t[1] extends past t[0]

	// check END
	int rete=3;
	f=0;
	s=1;
	if(tend[1]>tend[0]) { // f should be the longer transcript at the end
		f=1;
		s=0;
	}

	int j[2];
	j[0]=transfrag[t[0]]->nodes.Count()-1;
	j[1]=transfrag[t[1]]->nodes.Count()-1;

	len[2]=0;
	len[3]=0;
	while(transfrag[t[f]]->nodes[j[f]]>transfrag[t[s]]->nodes.Last()) { // s and f should have a node in common by now so fj doesn't have to go all the way back to 0
		if(no2gnode[transfrag[t[f]]->nodes[j[f]-1]]->end+1<no2gnode[transfrag[t[f]]->nodes[j[f]]]->start) { // found intron
			rete=1+f;
			if(rets<3 && rete!=rets) return 0;
			if(no2gnode[transfrag[t[s]]->nodes.Last()]->hardend) return 0; // there is no compatibility here
		}
		// to approximate:
		//len[2]+=no2gnode[transfrag[t[f]]->nodes[j[f]]]->len();
		// adj:
		if(j[f]<transfrag[t[f]]->nodes.Count()-1) len[2]+=no2gnode[transfrag[t[f]]->nodes[j[f]]]->len();
		else len[2]+=tend[f]-no2gnode[transfrag[t[f]]->nodes[j[f]]]->start+1;

		j[f]--;
	}
	if(transfrag[t[s]]->nodes[0]>transfrag[t[f]]->nodes[j[f]]) return 0; // no overlap (this shouldn't be the case since there is one node in common)

	if(j[f]==transfrag[t[f]]->nodes.Count()-1) len[2]+=tend[f]-tend[s];

	if(transfrag[t[s]]->nodes[j[s]]>transfrag[t[f]]->nodes[j[f]]) {
		if(!no2gnode[transfrag[t[f]]->nodes[j[f]]]->childpat[transfrag[t[s]]->nodes[j[s]]]) return 0; // I can not reach node j[s] from j[f]
		// len[3]=no2gnode[transfrag[t[s]]->nodes.Last()]->end-no2gnode[transfrag[t[f]]->nodes[j[f]]]->end; // to approximate dist
		//  adj if I want real distances:
		len[3]=tend[s]-no2gnode[transfrag[t[f]]->nodes[j[f]]]->end;
		while(transfrag[t[s]]->nodes[j[s]]>transfrag[t[f]]->nodes[j[f]]) {
			j[s]--;
			if(no2gnode[transfrag[t[s]]->nodes[j[s]]]->end+1<no2gnode[transfrag[t[s]]->nodes[j[s]+1]]->start) { // not continuous
				return 0;
			}
		}
	}
	else if(j[f]<transfrag[t[f]]->nodes.Count()-1) {
		len[2]+=no2gnode[transfrag[t[s]]->nodes.Last()]->end-tend[s];
	}

	if(i[f]>j[f] || i[s]>j[s]) return 0; // there should be i<=j

	if(f) len[2]=-len[2]; // a negative value signals that t[1] extends past t[0]

	// now if and is point at the same node in transcripts, and so do jf and js
	while(i[f]<=j[f] && i[s]<=j[s]) {
		if(transfrag[t[f]]->nodes[i[f]]==transfrag[t[s]]->nodes[i[s]]) { // skip equal nodes
			i[f]++;
			i[s]++;
		}
		else { // unequal
			if(transfrag[t[f]]->nodes[i[f]]>transfrag[t[s]]->nodes[i[s]]) { // make sure node f is always smallest
				int k=f;
				f=s;
				s=k;
			}
			while(transfrag[t[f]]->nodes[i[f]]<transfrag[t[s]]->nodes[i[s]] &&
					no2gnode[transfrag[t[f]]->nodes[i[f]-1]]->end+1==no2gnode[transfrag[t[f]]->nodes[i[f]]]->start) { // advance smallest node until I get to a gap
				i[f]++;
			}
			if(no2gnode[transfrag[t[f]]->nodes[i[f]-1]]->end+1==no2gnode[transfrag[t[f]]->nodes[i[f]]]->start) return 0; // gap in s is filled by transfrag f
			// see if gaps are both hard
			int *pos=gpos[edge(transfrag[t[f]]->nodes[i[f]-1],transfrag[t[f]]->nodes[i[f]],gno)];
			if(pos && transfrag[t[f]]->pattern[*pos]) { // this is a hard edge in t[f]
				pos=gpos[edge(transfrag[t[s]]->nodes[i[s]-1],transfrag[t[s]]->nodes[i[s]],gno)];
				if(pos && transfrag[t[s]]->pattern[*pos]) return 0; // edge is hard in both f and s
			}
			if(transfrag[t[f]]->nodes[i[f]]>transfrag[t[s]]->nodes[i[s]]) { // make sure node f is always smallest
				int k=f;
				f=s;
				s=k;
			}
			while(transfrag[t[f]]->nodes[i[f]]<transfrag[t[s]]->nodes[i[s]]) {
				if(no2gnode[transfrag[t[f]]->nodes[i[f]]]->end+1<no2gnode[transfrag[t[f]]->nodes[i[f]+1]]->start) return 0;
				i[f]++;
			}
		}
	}

	if(rets==1 || rete==1) return(1); // t[0] extends with introns past t[1] at at least one side
	if(rets==2 || rete==2) return(2); // t[1] extends with introns past t[0] at at least one side
	return(3); // t[0] and t[1] have compatible starts/ends
}

bool assign_incomplete_trf_to_nodes_APPLY_UNISPG(int t,int n1, int n2,GPVec<CGraphnodeUnispg>& no2gnode) {
	// a better version of this function could take into account also the fragments inner length
	bool assigned=false;

	for(int i=n1+1;i<n2;i++) {
		if(no2gnode[n1]->childpat[i] && no2gnode[i]->childpat[n2]) { // i is a child of n1, and n2 is a child of i
			assigned = binary_insert_trf_to_node_APPLY_UNISPG(t,no2gnode[i]->trf,0,no2gnode[i]->trf.Count()-1) or assigned;
		}
	}

	/* this version was taking a very long time
	for(int i=0;i<no2gnode[n1]->child.Count();i++) { // for all children of n1
		int child=no2gnode[n1]->child[i];
		if(child!=n2 && no2gnode[child]->childpat[n2]) { // n2 is child of n1's child
			binary_insert_trf_to_node_APPLY_UNISPG(t,no2gnode[child]->trf,0,no2gnode[child]->trf.Count()-1);
			assign_incomplete_trf_to_nodes_APPLY_UNISPG(t,child,n2,no2gnode);
		}
	}
	*/

	return assigned;
}


// this function inspects if there are multiple ways between two nodes in the transfrag; finds the first node like this and distributes the abundance of the transfrag between the two nodes in proportion to the edges that leave the first node
bool trf_real_APPLY_UNISPG(int t,GPVec<CGraphnodeUnispg>& no2gnode,GPVec<CTransfrag>& transfrag, GIntHash<int> &gpos,int gno) {

	int i=0; // index of the current node (the one I am looking at) in transfrag
	int nt=transfrag[t]->nodes.Count()-1; // index of last node in transfrag
	int n=transfrag[t]->nodes[i]; // start at first node in transfrag
	float totalcov=0;
	while(i<nt) {
		if(n+1==transfrag[t]->nodes[i+1]) { // if next node in transfrag is adjacent -> there is only one way to continue on transfrag path
			n=transfrag[t]->nodes[i+1];
			i++;
		}
		else { // there might be multiple ways to reach next node from node
			int *pos=gpos[edge(n,transfrag[t]->nodes[i+1],gno)];
			if(pos && transfrag[t]->pattern[*pos]) { // if there is an edge between node and nextnode in transfrag -> I don't have discontinuity there
				n=transfrag[t]->nodes[i+1];
				i++;
			}
			else { // found a potential discontinuity in transfrag pattern --> try to see if there are different branches going on from here
				CGraphnodeUnispg *node=no2gnode[n];
				int nextnode=transfrag[t]->nodes[i+1];
				int ntrf=node->trf.Count(); // total number of transfrags in node
				int np=0;  // upper bound on number of paths I can find from node to nextnode
				totalcov=0; // total abundance of transfrags going from node to children of node

				// I check if there are different paths to reach nextnode from node
				for(int j=0;j<node->child.Count();j++) {
					if(node->child[j]>nextnode) break; // if child of node is after nextnode -> I don't care
					if(node->child[j]==nextnode || no2gnode[node->child[j]]->childpat[nextnode]) { // I found a node on the path to nextnode
						CPath p(n,node->child[j]);
						int *pos=gpos[edge(n,node->child[j],gno)];
						if(pos) for(int f=0;f<ntrf;f++) {
							if(transfrag[node->trf[f]]->pattern[*pos]) { // this is a transfrag that goes through node and child; ATTENTION: it might not reach nextnode!! or be compatible with nextnode
								p.abundance+=transfrag[node->trf[f]]->abundance;
								totalcov+=transfrag[node->trf[f]]->abundance;
							}
						}
						if(p.abundance) {
							transfrag[t]->path.Add(p);
							np++;
						}
					}
				}
				if(!totalcov || !np) // there is no way to continue
					return true;
				if(np==1) { // only one way to go forward
					n=transfrag[t]->path[0].contnode;
					if(n==nextnode) i++;
					transfrag[t]->path.Clear();
				}
				else break; // break at first discontinue point in transfrag -> might be many other ways to continue but this functions doesn't care
			}
		}
	}
	if(i==nt) return true;

	// this implementation only assumes one break in transfrag
	for(int j=0;j<transfrag[t]->path.Count();j++) transfrag[t]->path[j].abundance/=totalcov;

	return false;
}


bool binary_insert_trf_to_node_APPLY_UNISPG(int t, GVec<int>& trf,int first,int last) {

	while(first<=last) {
		int mid=(first + last) / 2;  // compute mid point.
		if(t>trf[mid]) first=mid+1;
		else if(t<trf[mid]) last=mid-1;
		else return false; // transcript already inserted
	}

	int pos=first;

	trf.idxInsert(pos,t);
	return true;

}