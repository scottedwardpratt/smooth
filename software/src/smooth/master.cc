#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

CSmoothMaster::CSmoothMaster(){
	unsigned int NObs;
	parmap=new CparameterMap;
	parmap->ReadParsFromFile("smooth_data/Options/emulator_options.txt");
	int ranseed=parmap->getI("RANDY_SEED",time(NULL));
	randy=new Crandy(ranseed);
	
	string logfilename=parmap->getS("SmoothEmulator_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
   }
   SmoothEmulator_TrainingFormat=parmap->getS("SmoothEmulator_TrainingFormat","SMOOTH");
   SmoothEmulator_TestingFormat=parmap->getS("SmoothEmulator_TestingFormat","SMOOTH");
   string filename;
   
   filename="smooth_data/Info/observable_info.txt";
   
   observableinfo=new CObservableInfo(filename);
   NObs=observableinfo->NObservables;
   
	FullModelRunsDirName=parmap->getS("Smooth_FullModelRunsDirName","smooth_data/FullModelRuns");
   FullModelTestingRunsDirName=parmap->getS("Smooth_FullModelTestingRunsDirName","smooth_data/FullModelTestingRuns");
   
	SurmiseTrainingPointsFileName=parmap->getS("SmoothEmulator_SurmiseTrainingPointsFilename","smooth_data/SurmiseTrainingPoints.txt");
	SurmiseTrainingObsFileName=parmap->getS("SmoothEmulator_SurmiseTrainingObsFilename","smooth_data/SurmiseTrainingObs.txt");
   
   SurmiseTestingPointsFileName=parmap->getS("SmoothEmulator_SurmiseTestingPointsFilename","smooth_data/SurmiseTestingingPoints.txt");
   SurmiseTestingObsFileName=parmap->getS("SmoothEmulator_SurmiseTestingObsFilename","smooth_data/SurmiseTestingingObs.txt");
   
	filename="smooth_data/Info/prior_info.txt";
	priorinfo=new CPriorInfo(filename);
	NPars=priorinfo->NModelPars;
	parmap->set("SmoothEmulator_NPars",NPars);
	parmap->set("Smooth_NPars",NPars);
	CTrainingInfo::smoothmaster=this;
   CTestingInfo::smoothmaster=this;
	traininginfo = new CTrainingInfo(observableinfo,priorinfo);
   testinginfo = new CTestingInfo(observableinfo,priorinfo);
	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
   CSmoothEmulator::randy=randy;
   emulator.resize(NObs);
   
   for(unsigned int iy=0;iy<NObs;iy++){
      emulator[iy]=new CSmoothEmulator(observableinfo->observable_name[iy]);
   }
   ReadTrainingInfo();
   
}

void CSmoothMaster::CalcAllSigmaALambda(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->CalcSigmaALambda();
	}	
}

void CSmoothMaster::TuneAllY(){
	FILE *fptr=fopen("smooth_data/output_stuff/sigmalambda.txt","w");
	double sigmaAbar=0.0,Lambdabar=0.0;
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->Tune();
		fprintf(fptr,"%10.3f %10.5f\n",emulator[iY]->SigmaA,emulator[iY]->LAMBDA);
		sigmaAbar+=emulator[iY]->SigmaA;
		Lambdabar+=emulator[iY]->LAMBDA;
	}
	sigmaAbar=sigmaAbar/double(observableinfo->NObservables);
	Lambdabar=Lambdabar/double(observableinfo->NObservables);
	fclose(fptr);
}

void CSmoothMaster::TuneAllY(double LambdaSet){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      emulator[iY]->Tune(LambdaSet);
	}
}
	
void CSmoothMaster::TuneY(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(string obsname,double LambdaSet){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->Tune(LambdaSet);
}

void CSmoothMaster::TuneY(unsigned int iY){
	emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(unsigned int iY,double LambdaSet){
	emulator[iY]->Tune(LambdaSet);
}

void CSmoothMaster::GetY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
	emulator[iY]->GetYAndUncertainty(modelpars->Theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetY(unsigned int iY,vector<double> &theta,double &Y,double &SigmaY_emulator){
	emulator[iY]->GetYAndUncertainty(theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->GetYAndUncertainty(modelpars->Theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetY(string obsname,vector<double> &theta,double &Y,double &SigmaY_emulator){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->GetYAndUncertainty(theta,Y,SigmaY_emulator);
}

void CSmoothMaster::GetAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		GetY(iY,modelpars,Y[iY],SigmaY_emulator[iY]);
	}
}

void CSmoothMaster::GetAllY(vector<double> &Theta,vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		emulator[iY]->GetYAndUncertainty(Theta,Y[iY],SigmaY_emulator[iY]);
	}
}

double CSmoothMaster::GetYOnly(string obsname,CModelParameters *modelpars){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	double SigmaY_emulator,Y;
	emulator[iY]->GetYAndUncertainty(modelpars->Theta,Y,SigmaY_emulator);
	return Y;

}

double CSmoothMaster::GetYOnly(unsigned int iY,CModelParameters *modelpars){
	double SigmaY_emulator,Y;
	emulator[iY]->GetYAndUncertainty(modelpars->Theta,Y,SigmaY_emulator);
	return Y;
}

double CSmoothMaster::GetYOnly(string obsname,vector<double> &Theta){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	double SigmaY_emulator,Y;
	emulator[iY]->GetYAndUncertainty(Theta,Y,SigmaY_emulator);
	return Y;
}

double CSmoothMaster::GetYOnly(unsigned int iY,vector<double> &Theta){
	double SigmaY_emulator,Y;
	emulator[iY]->GetYAndUncertainty(Theta,Y,SigmaY_emulator);
	return Y;
}

double CSmoothMaster::GetYOnly(int iY,vector<double> Theta){
	double SigmaY_emulator,Y;
	emulator[iY]->GetYAndUncertainty(Theta,Y,SigmaY_emulator);
	return Y;
}

double CSmoothMaster::GetYOnlyPython(int DiY,vector<double> Theta){
	double SigmaY_emulator,Y;
	unsigned int iY=DiY;
	if(iY>=0 && iY<observableinfo->NObservables){
		emulator[iY]->GetYAndUncertainty(Theta,Y,SigmaY_emulator);
		return Y;
	}
	else
		return 0.0;
}

double CSmoothMaster::GetUncertainty(string obsname,vector<double> &Theta){
	double Y,Uncertainty;
	unsigned int iY=observableinfo->GetIPosition(obsname);
	emulator[iY]->GetYAndUncertainty(Theta,Y,Uncertainty);
	return Uncertainty;
}

double CSmoothMaster::GetUncertainty(unsigned int iY,vector<double> &Theta){
	double Y,Uncertainty;
	emulator[iY]->GetYAndUncertainty(Theta,Y,Uncertainty);
	return Uncertainty;
}

double CSmoothMaster::GetUncertainty(int iY,vector<double> Theta){
	double Y,Uncertainty;
	emulator[iY]->GetYAndUncertainty(Theta,Y,Uncertainty);
	return Uncertainty;
}

void CSmoothMaster::GetAllYOnly(CModelParameters *modelpars,vector<double> &Yvec){
	unsigned int NObservables=observableinfo->NObservables;
	Yvec.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		double Y,Uncertainty;
		emulator[iY]->GetYAndUncertainty(modelpars->Theta,Y,Uncertainty);
		Yvec[iY]=Y;
	}
}

void CSmoothMaster::GetAllYOnly(vector<double> &Theta,vector<double> &Yvec){
	unsigned int NObservables=observableinfo->NObservables;
	Yvec.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		double Y,Uncertainty;
		emulator[iY]->GetYAndUncertainty(Theta,Y,Uncertainty);
		Yvec[iY]=Y;
	}
}

vector<double> CSmoothMaster::GetYSigmaPython(int DiY,vector<double> theta){
	unsigned int iY=DiY;
	double Y,SigmaY_emulator;
	if(iY>=0 && iY<observableinfo->NObservables)
		emulator[iY]->GetYAndUncertainty(theta,Y,SigmaY_emulator);
	else{
		Y=SigmaY_emulator=0.0;
	}
	vector<double> YSigma;
	YSigma.resize(2);
	YSigma[0]=Y;
	YSigma[1]=SigmaY_emulator;
	return YSigma;		
}



