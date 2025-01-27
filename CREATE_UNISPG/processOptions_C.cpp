#include "processOptions_C.h"

// char* sprintTime() {
// 	static char sbuf[32];
// 	time_t ltime; /* calendar time */
// 	ltime=time(NULL);
// 	struct tm *t=localtime(&ltime);
// 	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
// 			t->tm_hour, t->tm_min, t->tm_sec);
// 	return(sbuf);
// }

void processCreateOptions(GArgs& args) {
	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n",VERSION);
	   exit(0);
	}

	if (args.getOpt("conservative")) {
	  isofrac=0.05;
	  singlethr=4.75;
	  readthr=1.5;
	  trim=false;
	}

	if (args.getOpt('t')) {
		trim=false;
	}

	if (args.getOpt("fr")) {
		fr_strand=true;
	}
	if (args.getOpt("rf")) {
		rf_strand=true;
		if(fr_strand) GError("Error: --fr and --rf options are incompatible.\n");
	}

	debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
	verbose=(args.getOpt('v')!=NULL);
	if (verbose) {
		fprintf(stderr, "Running StringTie " VERSION ". Command line:\n");
		args.printCmdLine(stderr);
	}

	keepTempFiles=(args.getOpt("keeptmp")!=NULL);

	//adaptive=!(args.getOpt('d')!=NULL);

	//complete=!(args.getOpt('i')!=NULL);
	includesource=!(args.getOpt('z')!=NULL);
	//EM=(args.getOpt('y')!=NULL);
	//weight=(args.getOpt('w')!=NULL);

	GStr s=args.getOpt('m');
	if (!s.is_empty()) {
	mintranscriptlen=s.asInt();
	if (mintranscriptlen<30)
		GError("Error: invalid -m value, must be >=30)\n");
	else if (mintranscriptlen<0) GError("Error: invalid -m value, must be >=0)\n");
	}

	multiMode=(args.getOpt("multi")!=NULL);
	if(multiMode) {
	} 

	graph_bed=(args.getOpt("graph_bed")!=NULL);

	s=args.getOpt("rseq");
	if (s.is_empty())
		s=args.getOpt('S');
	if (!s.is_empty()) {
		gfasta=new GFastaDb(s.chars());
	}

	//-- cram ref sequence
	s=args.getOpt("ref");
	if (s.is_empty())
		s=args.getOpt("cram-ref");
	if (!s.is_empty()) {
		cram_ref=s;
	}

	s=args.getOpt('x');
	if (!s.is_empty()) {
		//split by comma and populate excludeGSeqs
		s.startTokenize(" ,\t");
		GStr chrname;
		while (s.nextToken(chrname)) {
			excludeGseqs.Add(chrname.chars());
		}
	}

     /*
	 s=args.getOpt('n');
	 if (!s.is_empty()) {
		 sensitivitylevel=s.asInt();
		 if(sensitivitylevel<0) {
			 sensitivitylevel=0;
			 GMessage("sensitivity level out of range: setting sensitivity level at 0\n");
		 }
		 if(s>3) {
			 sensitivitylevel=3;
			 GMessage("sensitivity level out of range: setting sensitivity level at 2\n");
		 }
	 }
	*/


	s=args.getOpt('p');
	if (!s.is_empty()) {
		num_cpus=s.asInt();
		if (num_cpus<=0) num_cpus=1;
	}

	s=args.getOpt('a');
	if (!s.is_empty()) {
		junctionsupport=(uint)s.asInt();
	}

	s=args.getOpt('j');
	if (!s.is_empty()) junctionthr=s.asInt();

	s=args.getOpt('E');
	if (!s.is_empty()) sserror=s.asInt();

	s=args.getOpt('c');
	if (!s.is_empty()) {
		readthr=(float)s.asDouble();
		if (readthr<0.001) {
			GError("Error: invalid -c value, must be >=0.001)\n");
		}
	}


	s=args.getOpt('g');
	if (!s.is_empty()) {
		bundledist=s.asInt();
		if(bundledist>runoffdist) runoffdist=bundledist;
	}

	s=args.getOpt('F');
	if (!s.is_empty()) {
		fpkm_thr=(float)s.asDouble();
	}

	s=args.getOpt('T');
	if (!s.is_empty()) {
		tpm_thr=(float)s.asDouble();
	}

	s=args.getOpt('l');
	if (!s.is_empty()) label=s;

	s=args.getOpt('f');
	if (!s.is_empty()) {
		isofrac=(float)s.asDouble();
		if(isofrac>=1) GError("Miminum isoform fraction (-f coefficient: %f) needs to be less than 1\n",isofrac);
	}
	s=args.getOpt('M');
	if (!s.is_empty()) {
		mcov=(float)s.asDouble();
	}

	genefname=args.getOpt('A');
	if(!genefname.is_empty()) {
		geneabundance=true;
	}

	tmpfname=args.getOpt('o');

	// coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
	s=args.getOpt('s');
	if (!s.is_empty()) {
		singlethr=(float)s.asDouble();
		if (readthr<0.001) {
			GError("Error: invalid -s value, must be >=0.001)\n");
		}
	}

	if (args.getOpt('G')) {
	guidegff=args.getOpt('G');
	if (fileExists(guidegff.chars())>1)
		guided=true;
	else GError("Error: reference annotation file (%s) not found.\n",
				guidegff.chars());
	}
	s=args.getOpt("ptf");
	if (!s.is_empty()) {
	ptff=s;
	if (fileExists(ptff.chars())<=1)
		GError("Error: point features data file (%s) not found.\n",
				ptff.chars());
	}

	//  enableNames=(args.getOpt('E')!=NULL);

	retained_intron=(args.getOpt('i')!=NULL);

	//  isunitig=(args.getOpt('U')!=NULL);

	eonly=(args.getOpt('e')!=NULL);
	if(eonly) {
		eonly=false;
		includecov=true;
	}
	else if(eonly && !guided)
		GError("Error: invalid -e usage, GFF reference not given (-G option required).\n");

	 /* s=args->getOpt('P');
	 if (!s.is_empty()) {
		 if(!guided) GError("Error: option -G with reference annotation file has to be specified.\n");
		 c_out=fopen(s.chars(), "w");
		 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 partialcov=true;
	 }
	 else { */
		s=args.getOpt('C');
		if (!s.is_empty()) {
			if(!guided) GError("Error: invalid -C usage, GFF reference not given (-G option required).\n");
			c_out=fopen(s.chars(), "w");
			if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		}
	 //}
	int numbam=args.startNonOpt();
#ifndef GFF_DEBUG
	if (numbam < 1 ) {
	 	 GMessage("%s\nError: no input file provided!\n",USAGE);
	 	 exit(1);
	}
	if (args.startNonOpt()==0) {
        GMessage(USAGE);
        GMessage("\nError: no input file provided!\n");
        exit(1);
    }
#endif
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
	}
	else { // stdout
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

	
#ifdef B_DEBUG
	 GStr dbgfname(tmpfname);
	 dbgfname+=".dbg";
	 dbg_out=fopen(dbgfname.chars(), "w");
	 if (dbg_out==NULL) GError("Error creating debug output file %s\n", dbgfname.chars());
#endif

	tmpfname+=".tmp";
	f_out=fopen(tmpfname.chars(), "w");
	if (f_out==NULL) GError("Error creating output file %s\n", tmpfname.chars());

	nomulti=(args.getOpt('u')!=NULL);
}

