#include "msu_smooth/master.h"
using namespace std;

CSmoothMaster::CSmoothMaster(CparameterMap *parmap_set){
	parmap=parmap_set;
	int ranseed=parmap->getI("RANDY_SEED",-time(NULL));
	randy=new Crandy(ranseed);
	
	string filename=parmap->getS("OBSERVABLE_INFO_FILENAME","Info/observable_info.txt");
	observableinfo=new CObservableInfo(filename);
	
	ModelRunDirName=parmap->getS("SmoothEmulator_ModelRunDirName","model_runs");
	CoefficientsDirName=parmap->getS("SmoothEmulator_CoefficientsDirName","coefficients");
	
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

void CSmoothMaster::ReadTrainingInfo(){
	traininginfo->ReadTrainingInfo(ModelRunDirName);
	SetThetaTrain();
}

void CSmoothMaster::TuneAllY(){
	int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->Tune();
	}
}

void CSmoothMaster::TuneY(string obsname){
	int iY=observableinfo->GetIPosition(obsname);;
	emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(int iY){
	emulator[iY]->Tune();
}


void CSmoothMaster::GenerateCoefficientSamples(){
	int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		CLog::Info("Tuning Emulator for "+observableinfo->GetName(iY)+"\n");
		emulator[iY]->GenerateASamples();
	}
}

void CSmoothMaster::SetThetaTrain(){
	int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->SetThetaTrain();
	}
}

void CSmoothMaster::CalcY(int iY,CModelParameters *modelpars,double &Y,double &SigmaY){
	emulator[iY]->CalcY(modelpars,Y,SigmaY);
}

void CSmoothMaster::CalcY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY){
	int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcY(modelpars,Y,SigmaY);
}

void CSmoothMaster::CalcAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY){
	int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	SigmaY.resize(NObservables);
	for(int iY=0;iY<NObservables;iY++){
		CalcY(iY,modelpars,Y[iY],SigmaY[iY]);
	}
}

void CSmoothMaster::TestAtTrainingPts(){
	char pchars[CLog::CHARLENGTH];
	int itrain,iY;
	int NObservables=observableinfo->NObservables;
	vector<double> Y,SigmaY;
	Y.resize(NObservables);
	SigmaY.resize(NObservables);
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		for(iY=0;iY<NObservables;iY++){
			CalcY(iY,traininginfo->modelpars[itrain],Y[iY],SigmaY[iY]);
			snprintf(pchars,CLog::CHARLENGTH,
			"Y[%d]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,Y[iY],traininginfo->YTrain[iY][itrain],SigmaY[iY]);
			CLog::Info(pchars);
		}
	}
}
void CSmoothMaster::TestAtTrainingPts(int iY){
	char pchars[CLog::CHARLENGTH];
	int itrain;
	double Y,SigmaY;
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%d]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,Y,traininginfo->YTrain[iY][itrain],SigmaY);
		CLog::Info(pchars);
	}
}
void CSmoothMaster::TestAtTrainingPts(string obsname){
	char pchars[CLog::CHARLENGTH];
	int itrain,iY;
	double Y,SigmaY;
	iY=observableinfo->GetIPosition(obsname);
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%d]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,Y,traininginfo->YTrain[iY][itrain],SigmaY);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::WriteCoefficientsAllY(){
	for(int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->WriteCoefficients();
	}
}
void CSmoothMaster::WriteCoefficients(string obsname){
	int iY=observableinfo->GetIPosition(obsname);
	WriteCoefficients(iY);
}
void CSmoothMaster::WriteCoefficients(int iY){
	emulator[iY]->WriteCoefficients();
}

void CSmoothMaster::ReadCoefficientsAllY(){
	for(int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->ReadCoefficients();
	}
}
void CSmoothMaster::ReadCoefficients(string obsname){
	int iY=observableinfo->GetIPosition(obsname);
	ReadCoefficients(iY);
}
void CSmoothMaster::ReadCoefficients(int iY){
	emulator[iY]->ReadCoefficients();
}