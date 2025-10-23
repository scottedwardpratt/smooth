#include "msu_smooth/traininginfo.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;
CSmoothMaster* CTrainingInfo::smoothmaster=NULL;

CTrainingInfo::CTrainingInfo(CObservableInfo *observableinfo_set,CPriorInfo *priorinfo_set){
	observableinfo=observableinfo_set;
	priorinfo=priorinfo_set;
	CModelParameters::priorinfo=priorinfo;
	NObservables=observableinfo->NObservables;
	if(smoothmaster->SmoothEmulator_TrainingFormat == "SMOOTH"){
		ReadTrainingInfoSmoothFormat();
	}
	else if(smoothmaster->SmoothEmulator_TrainingFormat == "SURMISE"){
		string TrainingInfoFileName=smoothmaster->parmap->getS("SmoothEmulator_TrainingInfoFileName","traininginfo.txt");
		ReadTrainingInfoSurmiseFormat();
	}
	else{
		CLog::Fatal("SmoothEmulator_TrainingFormat not recognized,\n should be SMOOTH or SURMISE\n");
	}
	unsigned int iy,ntrain;
	YTrain.resize(NObservables);
	for(iy=0;iy<NObservables;iy++){
		YTrain[iy].resize(NTrainingPts);
		for(ntrain=0;ntrain<NTrainingPts;ntrain++)
			YTrain[iy][ntrain]=0.0;
	}

	CSmoothEmulator::NTrainingPts=NTrainingPts;

}

void CTrainingInfo::ReadTrainingInfoSmoothFormat(){
	unsigned int itrain,ilist,ifile,iy,nsuccess=0,ipar,nread;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char filename[300],obs_charname[300],mod_par_name[300];
	string obs_name;
	double y,x;
	FILE *fptr;
	
	if(smoothmaster->SmoothEmulator_TrainingFormat != "SMOOTH"){
		CLog::Fatal("SmoothEmulator_TrainingFormat should be set to SMOOTH\n if ReadTrainingInfo() is to be used\n");
	}
	//
	string NTrainingStr = smoothmaster->parmap->getS("SmoothEmulator_TrainingPts","1");
	vector<unsigned int> NTrainingList;
	NTrainingList.clear();
	stringstream ss(NTrainingStr);
	string token;
	string runfilename;
	int irun;
	bool exists;

	if(NTrainingStr=="all" || NTrainingStr=="All" || NTrainingStr=="ALL"){
		irun=0;
		exists=false;
		do{
			runfilename=smoothmaster->FullModelRunDirName+"/run"+to_string(irun);
			if(filesystem::exists(runfilename)){
				NTrainingList.push_back(irun);
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
					NTrainingList.push_back(i);
			}
			else {
				NTrainingList.push_back(stoi(token));
			}
		}
	}
	NTrainingPts = NTrainingList.size();
	
	modelpars.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]=new CModelParameters();
	}
	
	YTrain.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		YTrain[iy].resize(NTrainingPts);
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ifile=NTrainingList[itrain];
		snprintf(filename,300,"%s/run%u/obs.txt",smoothmaster->FullModelRunDirName.c_str(),ifile);
		fptr=fopen(filename,"r");
		nsuccess=0;
		do{
			fscanf(fptr,"%s",obs_charname);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&y);
				obs_name=string(obs_charname);
				iy=smoothmaster->observableinfo->GetIPosition(obs_name);
				YTrain[iy][itrain]=y;
				nsuccess+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}

	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTrainingInfo::ReadTrainInfo, only read in "+to_string(nsuccess)+" observables from file "+string(filename)+"\n");

	for(itrain=0;itrain<NTrainingPts;itrain++){
		ilist=NTrainingList[itrain];
		snprintf(filename,300,"%s/run%u/model_parameters.txt",smoothmaster->FullModelRunDirName.c_str(),ilist);
		fptr=fopen(filename,"r");
		nread=0;
		do{
			fscanf(fptr,"%s",mod_par_name);
			if(!feof(fptr)){
				fscanf(fptr,"%lf",&x);
				ipar=priorinfo->GetIPosition(mod_par_name);
				modelpars[itrain]->X[ipar]=x;
				nread+=1;
			}
		}while(!feof(fptr));
		fclose(fptr);
	}
	if(nread!=priorinfo->NModelPars){
		CLog::Fatal("Only read in "+to_string(nread)+" parameter values from "+string(filename)+". But there are "+to_string(priorinfo->NModelPars)+" parameters needed.\n");
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
	}

}

void CTrainingInfo::ReadTrainingInfoSurmiseFormat(){
	if(smoothmaster->SmoothEmulator_TrainingFormat != "SURMISE"){
		CLog::Fatal("SmoothEmulator_TrainingFormat should be set to SURMISE\n if ReadTrainingInfoSurmiseFormat() is to be used\n");
	}
	unsigned int itrain,ipar,iobs;
	unsigned int NModelPars=CModelParameters::NModelPars;
	unsigned int NObs=smoothmaster->observableinfo->NObservables;
	char dummy[10000];
	string obs_name,filename;
	double y,x;
	
	filename=smoothmaster->SurmiseTrainingPointsFileName;
	FILE *fptr=fopen(filename.c_str(),"r");
	itrain=0;
	do{
		for(ipar=0;ipar<NModelPars;ipar++){
			fscanf(fptr,"%lf",&x);
			if(!feof(fptr)){
				if(ipar==0){
					modelpars.push_back(NULL);
					modelpars[itrain]=new CModelParameters();
				}
				modelpars[itrain]->X[ipar]=x;
			}
		}
		fgets(dummy,10000,fptr);
		if(!feof(fptr))
			itrain+=1;
	}while(!feof(fptr));
	fclose(fptr);
	
	NTrainingPts=itrain;
	filename=smoothmaster->SurmiseTrainingObsFileName;
	fptr=fopen(filename.c_str(),"r");
	
	YTrain.resize(NObs);
	for(iobs=0;iobs<NObs;iobs++){
		YTrain[iobs].resize(NTrainingPts);
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(iobs=0;iobs<NObs;iobs++){
			fscanf(fptr,"%lf",&y);
			if(feof(fptr)){
				CLog::Fatal("reading training info: not enough lines in "+smoothmaster->SurmiseTrainingObsFileName+"\n");
			}
			YTrain[iobs][itrain]=y;
		}
		fgets(dummy,10000,fptr);
	}
	
	fclose(fptr);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelpars[itrain]->TranslateX_to_Theta();
	}
	
}

