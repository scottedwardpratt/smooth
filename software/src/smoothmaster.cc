#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
using namespace std;

CSmoothMaster::CSmoothMaster(CparameterMap *parmap_set){
	parmap=parmap_set;
	int ranseed=parmap->getI("RANDY_SEED",-time(NULL));
	randy=new Crandy(ranseed);
	
	string filename=parmap->getS("OBSERVABLE_INFO_FILENAME","observable_info.txt");
	observableinfo=new CObservableInfo(filename);
	
	filename=parmap->getS("PRIOR_INFO_FILENAME","prior_info.txt");
	priorinfo=new CPriorInfo(filename);
	
	traininginfo=new CTrainingInfo(observablinfo,priorinfo);
	
	emulator.resize(observableinfo->NObservables);
	
	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	CTrainingInfo::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
	CSmoothEmulator::randy=randy;
	
	
	for(int i=0;i<observableinfo->NObservables;i++){
		emulator[i]=new CSmoothEmulator(observableinfo->observable_name[i]);
	}
	
	
	smooth=new CSmooth(parmap);
	
	NPars=parmap->getD("SmoothEmulator_NPars",0);
	parmap->set("Smooth_NPars",NPars);
	
	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	
	
	
}
