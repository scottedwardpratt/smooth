#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

CSmoothMaster::CSmoothMaster(CparameterMap *parmap_set){
	parmap=parmap_set;
	int ranseed=parmap->getI("RANDY_SEED",time(NULL));
	randy=new Crandy(ranseed);
	
	printf("howdy a\n");
	
	string logfilename=parmap->getS("SmoothEmulator_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	
	printf("howdy b\n");

	string filename;
	UsePCA=parmap->getB("SmoothEmulator_UsePCA",false);
	if(UsePCA){
		filename="Info/pca_info.txt";
		CoefficientsDirName="coefficients_pca";
	}
	else{
		filename="Info/observable_info.txt";
		CoefficientsDirName="coefficients";
	}
	printf("howdy bb, filename=%s\n",filename.c_str());
	observableinfo=new CObservableInfo(filename);
	
	printf("howdy bbb\n");

	ModelRunDirName=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	
	printf("howdy c\n");

	filename="Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(filename);
	NPars=priorinfo->NModelPars;
	parmap->set("SmoothEmulator_NPars",NPars);
	parmap->set("Smooth_NPars",NPars);

	string NTrainingStr = parmap->getS("SmoothEmulator_TrainingPts","1");
	
	printf("howdy f\n");
	
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
	CSmoothEmulator::smooth=smooth;
	CSmoothEmulator::smoothmaster=this;
	emulator.resize(observableinfo->NObservables);
	
	printf("howdy g\n");
	
	for(unsigned int i=0;i<observableinfo->NObservables;i++){
		emulator[i]=new CSmoothEmulator(observableinfo->observable_name[i]);
	}

}

void CSmoothMaster::TuneAllY(){
	printf("check aa\n");
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->Tune();
	}
}

void CSmoothMaster::TuneY(string obsname){
	unsigned int iY=observableinfo->GetIPosition(obsname);
	if(iY==3)
		emulator[3]->A[0]=5.5;
	emulator[iY]->Tune();
}

void CSmoothMaster::TuneY(unsigned int iY){
	if(iY==3)
		emulator[3]->A[0]=5.5;
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

void CSmoothMaster::CalcYdYdTheta(unsigned int iY,CModelParameters *modelpars,double &Y,
double &SigmaY_emulator,vector<double> &dYdTheta){
	emulator[iY]->CalcYDYDTheta(modelpars,Y,dYdTheta,SigmaY_emulator);
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

void CSmoothMaster::CalcAllYdYdTheta(CModelParameters *modelpars,vector<double> &Y,
vector<double> &SigmaY_emulator,vector<vector<double>> &dYdTheta){
	unsigned int NObservables=observableinfo->NObservables;
	Y.resize(NObservables);
	dYdTheta.resize(NObservables);
	SigmaY_emulator.resize(NObservables);
	for(unsigned int iY=0;iY<NObservables;iY++){
		CalcYdYdTheta(iY,modelpars,Y[iY],SigmaY_emulator[iY],dYdTheta[iY]);
		dYdTheta.resize(NPars);
	}
}

void CSmoothMaster::CalcAllLogP(){
	unsigned int NObservables=observableinfo->NObservables;
	double logP,logPbar=0.0,A2barRatioBar=0.0,SigmaAbar=0.0;
	for(unsigned int iY=0;iY<NObservables;iY++){
		emulator[iY]->CalcExactLogP();
		logP=emulator[iY]->logP;
		CLog::Info("iY="+to_string(iY)+": logP="+to_string(logP)+", A2Ratio="+to_string(emulator[iY]->A2barRatio)+"\n");
		logPbar+=logP;
		A2barRatioBar+=emulator[iY]->A2barRatio;
		SigmaAbar+=emulator[iY]->SigmaA;
	}
	logPbar=logPbar/double(NObservables);
	A2barRatioBar=A2barRatioBar/double(NObservables);
	SigmaAbar=SigmaAbar/double(NObservables);
	//CLog::Info("logPbar="+to_string(logPbar)+", A2barRatio="+to_string(A2barRatioBar)+"\n");
	CLog::Info(to_string(emulator[0]->LAMBDA)+" "+to_string(logPbar)+" "+to_string(A2barRatioBar)
		+" "+to_string(SigmaAbar)+"\n");
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

void CSmoothMaster::TestVsFakeModel(){
	char pchars[CLog::CHARLENGTH];
	unsigned int iY,ifake,nfake=100,ipar;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator,fakeY;
	CModelParameters fakepars[nfake];
	CLog::Info("--- TESTING VS FAKE MODEL ----\n");
	FILE *fptr,*fptr_out;
	string filename;
	
	
	
	for(iY=0;iY<NObservables;iY++){
		printf("SigmaA[%d]=%g\n",iY,emulator[iY]->SigmaA);
		filename="fakedata/"+observableinfo->observable_name[iY]+".txt";
		fptr=fopen(filename.c_str(),"r");
		filename="fakedata/fake_"+observableinfo->observable_name[iY]+".txt";
		fptr_out=fopen(filename.c_str(),"w");
		
		for(ifake=0;ifake<nfake;ifake++){
			for(ipar=0;ipar<NPars;ipar++){
				fscanf(fptr,"%lf",&fakepars[ifake].Theta[ipar]);
			}
			fscanf(fptr,"%lf",&fakeY);
			CalcY(iY,&fakepars[ifake],Y,SigmaY_emulator);
			snprintf(pchars,CLog::CHARLENGTH,
			"Y[%u]=%10.3e =? %10.3e,    SigmaY_emulator=%12.5e\n",
			iY,Y,fakeY,SigmaY_emulator);
			fprintf(fptr_out,"%12.5e  %12.5e %12.5e\n",fakeY,Y,SigmaY_emulator);	
		}
		fclose(fptr);		
		fclose(fptr_out);
	}
}

