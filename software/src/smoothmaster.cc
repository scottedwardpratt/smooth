#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
using namespace std;

CSmoothMaster::CSmoothMaster(CparameterMap *parmap_set){
	parmap=parmap_set;
	string obs_info_filename=parmap->getS("SMOOTH_OBSERVABLE_INFO_FILENAME","observable_info.txt");
	observableinfo=new CObservableInfo(obs_info_filename);
	emulator.resize(observableinfo->NObservables);
	
	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
	
	
	for(int i=0;i<observableinfo->NObservables;i++){
		emulator[i]=new CSmoothEmulator(observableinfo);
	}
	int ranseed=parmap->getI("RANDY_SEED",-time(NULL));
	randy=new Crandy(ranseed);
	smooth=new CSmooth(parmap);
	
	NPars=parmap->getD("SmoothEmulator_NPars",0);
	parmap->set("Smooth_NPars",NPars);
	
	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	
	
	
}
