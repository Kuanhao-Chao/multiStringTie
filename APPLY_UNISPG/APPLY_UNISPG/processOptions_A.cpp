#include "processOptions_A.h"

void processApplyOptions(GArgs& args) {

    fprintf(stderr, "Inside 'processApplyOptions'\n");
	debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n",VERSION);
	   exit(0);
	}	
	GStr s=args.getOpt("dot");
	if (s.is_empty()) {
	   s=args.getOpt("d");
	   fprintf(stderr, "It's empty%s \n", s.chars());
	}
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
	
	tmpfname=args.getOpt('o');

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
