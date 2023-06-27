#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
using namespace std;

CSmoothMaster::CSmoothMaster(CparameterMap *parmap_set){
	parmap=parmap_set;
	int ranseed=parmap->getI("RANDY_SEED",-time(NULL));
	randy=new Crandy(ranseed);
	
	string filename=parmap->getS("OBSERVABLE_INFO_FILENAME","Info/observable_info.txt");
	observableinfo=new CObservableInfo(filename);
	
	filename=parmap->getS("PRIOR_INFO_FILENAME","Info/prior_info.txt");
	priorinfo=new CPriorInfo(filename);
	NPars=priorinfo->NModelPars;
	parmap->set("Smooth_NPars",NPars);
	
	NTrainingPts=parmap->getI("SmoothEmulator_NTrainingPts",0);
	CTrainingInfo::smoothmaster=this;
	traininginfo=new CTrainingInfo(NTrainingPts,observableinfo,priorinfo);
	
	smooth=new CSmooth(parmap);

	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
	CSmoothEmulator::randy=randy;
	CSmoothEmulator::NTrainingPts=NTrainingPts;
	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	emulator.resize(observableinfo->NObservables);
	for(int i=0;i<observableinfo->NObservables;i++){
		emulator[i]=new CSmoothEmulator(observableinfo->observable_name[i]);
	}
	
}

void CSmoothMaster::ReadTrainingInfo(string modeldir){
	traininginfo->ReadTrainingInfo(modeldir);
	SetThetaTrain();
}

void CSmoothMaster::TuneA(){
	int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->TuneA();
	}
}

void CSmoothMaster::GenerateASamples(){
	int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->GenerateASamples();
	}
}

void CSmoothMaster::SetThetaTrain(){
	int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->SetThetaTrain();
	}
}