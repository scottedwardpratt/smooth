#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;


CSmoothMaster::CSmoothMaster(CparameterMap *parmap_set){
	parmap=parmap_set;
	int ranseed=parmap->getI("RANDY_SEED",time(NULL));
	randy=new Crandy(ranseed);
	
	string logfilename=parmap->getS("SmoothEmulator_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}

	string filename;
	if(UsePCA){
		filename="Info/pca_info.txt";
		CoefficientsDirName="coefficients_pca";
	}
	else{
		filename="Info/observable_info.txt";
		CoefficientsDirName="coefficients";
	}
	observableinfo=new CObservableInfo(filename);

	ModelRunDirName=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");

	filename="Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(filename);
	NPars=priorinfo->NModelPars;
	parmap->set("SmoothEmulator_NPars",NPars);
	parmap->set("Smooth_NPars",NPars);
	UsePCA=parmap->getB("SmoothEmulator_UsePCA",false);

	string NTrainingStr = parmap->getS("SmoothEmulator_TrainingPts","1");
	
	vector<unsigned int> NTrainingList;
	stringstream ss(NTrainingStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos) {

			unsigned int start = stoi(token.substr(0, pos));
			unsigned int end = stoi(token.substr(pos+1));

			for (unsigned int i = start; i <= end; i++)
			NTrainingList.push_back(i);
		}
		else {

			NTrainingList.push_back(stoi(token));
		}
	}

	CTrainingInfo::smoothmaster=this;
	traininginfo = new CTrainingInfo(NTrainingList,observableinfo,priorinfo);

	smooth=new CSmooth(parmap);

	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
	CSmoothEmulator::randy=randy;
	CSmoothEmulator::NTrainingPts=NTrainingList.size();
	printf("NTrainingPts=%u\n",CSmoothEmulator::NTrainingPts);
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	emulator.resize(observableinfo->NObservables);
	
	for(unsigned int i=0;i<observableinfo->NObservables;i++){
		emulator[i]=new CSmoothEmulator(observableinfo->observable_name[i]);
	}

}

void CSmoothMaster::ReadTrainingInfo(){
	traininginfo->ReadTrainingInfo(ModelRunDirName);
	SetThetaTrain();
}

void CSmoothMaster::TuneAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->Tune();
	}
}

void CSmoothMaster::TuneY(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);;
	emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(unsigned int iY){
	emulator[iY]->Tune();
}

void CSmoothMaster::GenerateCoefficientSamples(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		CLog::Info("Tuning Emulator for "+observableinfo->GetName(iY)+"\n");
		emulator[iY]->GenerateASamples();
	}
}

void CSmoothMaster::SetThetaTrain(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->SetThetaTrain();
	}
}

void CSmoothMaster::CalcY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
	emulator[iY]->CalcY(modelpars,Y,SigmaY_emulator);
}

void CSmoothMaster::CalcY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->CalcY(modelpars,Y,SigmaY_emulator);
}

void CSmoothMaster::CalcAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcY(iY,modelpars,Y[iY],SigmaY_emulator[iY]);
	}
}

void CSmoothMaster::TestAtTrainingPts(){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain,iY;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator;
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		for(iY=0;iY<NObservables;iY++){
			CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
			snprintf(pchars,CLog::CHARLENGTH,
			"Y[%u]=%10.3e =? %10.3e,    SigmaY_emulator=%12.5e\n",iY,Y,traininginfo->YTrain[iY][itrain],SigmaY_emulator);
			CLog::Info(pchars);
		}
	}
}

void CSmoothMaster::TestAtTrainingPts(unsigned int iY){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain;
	double Y,SigmaY_emulator;
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,Y,traininginfo->YTrain[iY][itrain],SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::TestAtTrainingPts(string obsname){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain,iY;
	double Y,SigmaY_emulator;
	iY=observableinfo->GetIPosition(obsname);
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		CalcY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,Y,traininginfo->YTrain[iY][itrain],SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::WriteCoefficientsAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->WriteCoefficients();
	}
}

void CSmoothMaster::WriteCoefficients(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	WriteCoefficients(iY);
}

void CSmoothMaster::WriteCoefficients(unsigned int iY){
	emulator[iY]->WriteCoefficients();
}

void CSmoothMaster::ReadCoefficientsAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->ReadCoefficients();
	}
}

void CSmoothMaster::ReadCoefficients(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	ReadCoefficients(iY);
}

void CSmoothMaster::ReadCoefficients(unsigned int iY){
	emulator[iY]->ReadCoefficients();
}
