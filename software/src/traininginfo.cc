#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

using namespace std;

CTrainingInfo::CTrainingInfo(int NTrainingPts_set,CObservableInfo *observableinfo,CPriorInfo *priorinfo){
	NObservables=observableinfo->NObservables;
	NTrainingPts=NTrainingPts_set;
	int iy,ntrain;
	YTrain.resize(NObservables);
	SigmaYTrain.resize(NObservables);
	for(iy=0;iy<NObservables;iy++){
		YTrain[iy].resize(NTrainingPts);
		SigmaYTrain[iy].resize(NTrainingPts);
		for(ntrain=0;ntrain<NTrainingPts;ntrain++)
			YTrain[iy][ntrain]=SigmaYTrain[iy][ntrain]=0.0;
	}

	modelpars.resize(NTrainingPts);
	for(ntrain=0;ntrain<NTrainingPts;ntrain++){
		modelpars[ntrain]=new CModelParameters(smoothmaster->priorinfo);
	}

}

void CTrainingInfo::ReadTrainInfo(string rundirname){
	bool success;
	int irun,iy,nsuccess=0;
	char filename[300],obs_charname[300];
	string obs_name;
	double y,sigmay;
	FILE *fptr;
	for(irun=0;irun<NTrainingPts;irun++){
		success=false;
		snprintf(filename,300,"%s/run%d/obs.txt",rundirname.c_str(),irun);
		fptr=fopen(filename,"r");
		do{
			fscanf(fptr,"%s %lf %lf",obs_charname,&y,&sigmay);
			obs_name=string(obs_charname);
			iy=smoothmaster->observableinfo->GetIPosition(obs_name);
			if(iy<0)
				CLog::Fatal("In CTrainingInfo::ReadYTrain, reading filename"+string(filename)+", but cannot recognize observable "+obs_name+"\n");
			YTrain[iy][irun]=y;
			SigmaYTrain[iy][irun]=sigmay;
			nsuccess+=1;			
		}while(!feof(fptr) && success==false);
		if(success==false){
			CLog::Fatal("For irun="+to_string(irun)+" cannot find YTrain for observable"+obs_name+"\n");
		}
	}
	if(nsuccess!=smoothmaster->observableinfo->NObservables)
		CLog::Fatal("In CTrainingInfo::ReadTrainInfo, only read in "+to_string(nsuccess)+"observables from file "+string(filename)+"\n");
}
