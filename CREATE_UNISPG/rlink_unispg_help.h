void parse_trf_unispg(int maxi,int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,bool first,int strand,GList<CPrediction>& pred,GVec<float>& nodecov,
		GBitVec& istranscript,GBitVec& usednode,float maxcov,GBitVec& prevpath);

int guides_pushmaxflow_unispg(int gno,int edgeno,GIntHash<int>& gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,GVec<CGuide>& guidetrf,int& geneno,
		int s,GList<CPrediction>& pred,GVec<float>& nodecov,GBitVec& istranscript,GBitVec& pathpat,bool &first,GPVec<GffObj>& guides,GVec<int> &guidepred, BundleData *bdata);

void get_trf_long_mix_unispg(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,int& geneno,int strand,
		GList<CPrediction>& pred,GVec<int>& trflong,GVec<float>& nodecovall,GBitVec& istranscript,GBitVec& prevpath,BundleData *bdata,bool &first);

void get_trf_long_unispg(int gno,int edgeno, GIntHash<int> &gpos,GPVec<CGraphnode>& no2gnode,GPVec<CTransfrag>& transfrag,
		int& geneno,int strand,GList<CPrediction>& pred,GVec<int>& trflong,BundleData *bdata);