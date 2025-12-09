#include "msu_smooth/testinginfo.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;
CSmoothMaster* CTestingInfo::smoothmaster=NULL;

CTestingInfo::CTestingInfo(CObservableInfo *observableinfo_set,CPriorInfo *priorinfo_set){
	observableinfo=observableinfo_set;
	priorinfo=priorinfo_set;
	CModelParameters::priorinfo=priorinfo;
	NObservables=observableinfo->NObservables;
	if(smoothmaster->SmoothEmulator_TestingFormat == "SMOOTH"){
		ReadTestingInfoSmoothFormat();
	}
	else if(smoothmaster->SmoothEmulator_TestingFormat == "SURMISE"){
		string TestingInfoFileName=smoothmaster->parmap->getS("SmoothEmulator_TestingInfoFileName","testinginfo.txt");
		ReadTestingInfoSurmiseFormat();
	}
	else{
		CLog::Fatal("SmoothEmulator_TestingFormat not recognized,\n should be SMOOTH or SURMISE\n");
	}
	unsigned int iy,ntest;
	YTest.resize(NObservables);
	for(iy=0;iy<NObservables;iy++){
		YTest[iy].resize(NTestingPts);
		for(ntest=0;ntest<NTestingPts;ntest++)
			YTest[iy][ntest]=0.0;
	}

	CSmoothEmulator::NTestingPts=NTestingPts;

}

void CTestingInfo::ReadTestingInfoSmoothFormat(){
	unsigned int itest,ilist,ifile,iy,nsuccess=0,ipar,nread;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char filename[300],obs_charname[300],mod_par_name[300];
	string obs_name;
	double y,x;
	FILE *fptr;	if(smoothmaster->SmoothEmulator_TestingFormat != "SMOOTH"){
		CLog::Fatal("SmoothEmulator_TestingFormat should be set to SMOOTH\n if ReadTestingInfo() is to be used\n");
	}
	//
	string NTestingStr = smoothmaster->parmap->getS("SmoothEmulator_TestingPts","all");
	vector<unsigned int> NTestingList;
	NTestingList.clear();
	stringstream ss(NTestingStr);
	string token;
	string runfilename;
   
	int irun;
	bool exists;
   
	if(NTestingStr=="all" || NTestingStr=="All" || NTestingStr=="ALL"){
		irun=0;
		exists=false;
		do{
			runfilename=smoothmaster->FullModelTestingRunsDirName+"/run"+to_string(irun);
			if(filesystem::exists(runfilename)){
				NTestingList.push_back(irun);
				exists=true;
			}
			else
				exists=false;
			irun+=1;
		}while(exists);
	}
	else{
		while(getline(ss, token, ',')) {
			size_t pos = token.find("-");
			if (pos != string::npos){

				unsigned int start = stoi(token.substr(0, pos));
				unsigned int end = stoi(token.substr(pos+1));

				for (unsigned int i = start; i <= end; i++)
					NTestingList.push_back(i);
			}
			else {
				NTestingList.push_back(stoi(token));
			}
		}
	}
	NTestingPts=NTestingList.size();
   	
	modelpars.resize(NTestingPts);
	for(itest=0;itest<NTestingPts;itest++){
		modelpars[itest]=new CModelParameters();
	}
	YTest.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		YTest[iy].resize(NTestingPts);
	}
   
	for(itest=0;itest<NTestingPts;itest++){
		ifile=NTestingList[itest];
		snprintf(filename,300,"%s/run%u/obs.txt",smoothmaster->FullModelTestingRunsDirName.c_str(),ifile);
		fptr=fopen(filename,"r");
		nsuccess=0;
		do{
			fscanf(fptr,"%s",obs_charname);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&y);
				obs_name=string(obs_charname);
				iy=smoothmaster->observableinfo->GetIPosition(obs_name);
				YTest[iy][itest]=y;
				nsuccess+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}

	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTestingInfo::ReadTestInfo, only read in "+to_string(nsuccess)+" observables from file "+string(filename)+"\n");

	for(itest=0;itest<NTestingPts;itest++){
		ilist=NTestingList[itest];
		snprintf(filename,300,"%s/run%u/model_parameters.txt",smoothmaster->FullModelTestingRunsDirName.c_str(),ilist);
		fptr=fopen(filename,"r");
		nread=0;
		do{
			fscanf(fptr,"%s",mod_par_name);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&x);
				ipar=priorinfo->GetIPosition(mod_par_name);
				modelpars[itest]->X[ipar]=x;
				nread+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}
	if(nread!=priorinfo->NModelPars){
		CLog::Fatal("Only read in "+to_string(nread)+" parameter values from "+string(filename)+". But there are "+to_string(priorinfo->NModelPars)+" parameters needed.\n");
	}
	
	for(itest=0;itest<NTestingPts;itest++){
		modelpars[itest]->TranslateX_to_Theta();
	}

}

void CTestingInfo::ReadTestingInfoSurmiseFormat(){
	if(smoothmaster->SmoothEmulator_TestingFormat != "SURMISE"){
		CLog::Fatal("SmoothEmulator_TestingFormat should be set to SURMISE\n if ReadTestingInfoSurmiseFormat() is to be used\n");
	}
	unsigned int itest,ipar,iobs;
	unsigned int NModelPars=CModelParameters::NModelPars;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char dummy[10000];
	string obs_name,filename;
	double y,x;
	
	filename=smoothmaster->SurmiseTestingParsFileName;
	FILE *fptr=fopen(filename.c_str(),"r");
	itest=0;
	do{
		for(ipar=0;ipar<NModelPars;ipar++){
			fscanf(fptr,"%lf",&x);
			if(!feof(fptr)){
				if(ipar==0){
					modelpars.push_back(NULL);
					modelpars[itest]=new CModelParameters();
				}
				modelpars[itest]->X[ipar]=x;
			}
		}
		fgets(dummy,10000,fptr);
		if(!feof(fptr))
			itest+=1;
	}while(!feof(fptr));
	fclose(fptr);
	
	NTestingPts=itest;
	filename=smoothmaster->SurmiseTestingObsFileName;
	fptr=fopen(filename.c_str(),"r");
	
	YTest.resize(NObs);
	for(iobs=0;iobs<NObs;iobs++){
		YTest[iobs].resize(NTestingPts);
	}
	
	for(itest=0;itest<NTestingPts;itest++){
		for(iobs=0;iobs<NObs;iobs++){
			fscanf(fptr,"%lf",&y);
			if(feof(fptr)){
				CLog::Fatal("reading testing info: not enough lines in "+smoothmaster->SurmiseTestingObsFileName+"\n");
			}
			YTest[iobs][itest]=y;
		}
		fgets(dummy,10000,fptr);
	}
	
	fclose(fptr);
	
	for(itest=0;itest<NTestingPts;itest++){
		modelpars[itest]->TranslateX_to_Theta();
	}
	
}
