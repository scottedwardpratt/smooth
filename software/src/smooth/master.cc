#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

CSmoothMaster::CSmoothMaster(){
	unsigned int NObs;
	unsigned int iZ;
	CPCA *pca=NULL;
	parmap=new CparameterMap;
	parmap->ReadParsFromFile("smooth_data/smooth_parameters/emulator_parameters.txt");
	int ranseed=parmap->getI("RANDY_SEED",time(NULL));
	randy=new Crandy(ranseed);
	
	string logfilename=parmap->getS("SmoothEmulator_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	SmoothEmulator_TrainingFormat=parmap->getS("SmoothEmulator_TrainingFormat","training_format_smooth");
	string filename;
	UsePCA=parmap->getB("SmoothEmulator_UsePCA",false);
	if(UsePCA){
		filename="smooth_data/PCA_Info/observable_info.txt";
		CoefficientsDirName="coefficients_pca";
		pca=new CPCA();
	}
	else{
		filename="smooth_data/Info/observable_info.txt";
		CoefficientsDirName="coefficients";
	}
	observableinfo=new CObservableInfo(filename);
	NObs=observableinfo->NObservables;
	
	ModelRunDirName=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	TrainingThetasFileName=parmap->getS("SmoothEmulator_TrainingThetasFilename","TrainingThetas.txt");
	TrainingObsFileName=parmap->getS("SmoothEmulator_TrainingObsFilename","TrainingObs.txt");
	
	filename="smooth_data/Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(filename);
	NPars=priorinfo->NModelPars;
	parmap->set("SmoothEmulator_NPars",NPars);
	parmap->set("Smooth_NPars",NPars);
	CTrainingInfo::smoothmaster=this;
	traininginfo = new CTrainingInfo(observableinfo,priorinfo);

	CSmoothEmulator::NPars=NPars;
	CSmoothEmulator::smoothmaster=this;
	CSmoothEmulator::parmap=parmap;
	CSmoothEmulator::randy=randy;
	emulator.resize(NObs);
	
	pca_ignore.resize(NObs);
	if(UsePCA){
		pca->ReadTransformationInfo();
		pca_minvariance=parmap->getD("SmoothEmulator_PCAMinVariance",0.0);
		for(iZ=0;iZ<NObs;iZ++)
			pca_ignore[iZ]=false;
		for(iZ=0;iZ<NObs;iZ++){
			if(pca->eigvals(iZ)<pca_minvariance)
				pca_ignore[iZ]=true;
			else
				pca_ignore[iZ]=false;
		}
		delete pca;
	}
	else{
		for(iZ=0;iZ<NObs;iZ++)
			pca_ignore[iZ]=false;
	}
	for(unsigned int iy=0;iy<NObs;iy++){
		if(!UsePCA || !pca_ignore[iy]){
			emulator[iy]=new CSmoothEmulator(observableinfo->observable_name[iy]);
		}
	}
	ReadTrainingInfo();
	
}

void CSmoothMaster::CalcAllSigmaALambda(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->CalcSigmaALambda();
	}	
}

void CSmoothMaster::TuneAllY(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		if((UsePCA && !pca_ignore[iY]) || !UsePCA){
			//CLog::Info("---- Tuning for "+observableinfo->observable_name[iY]+" ----\n");
			emulator[iY]->Tune();
		}
	}
}

void CSmoothMaster::TuneAllY(double LambdaSet){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		if((UsePCA && !pca_ignore[iY]) || !UsePCA){
			//CLog::Info("---- Tuning for "+observableinfo->observable_name[iY]+" ----\n");
			emulator[iY]->Tune(LambdaSet);
		}
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

/*void CSmoothMaster::GenerateCoefficientSamples(){
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		CLog::Info("Tuning Emulator for "+observableinfo->GetName(iY)+"\n");
		emulator[iY]->GenerateASamples();
	}
}*/

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

void CSmoothMaster::TestAtTrainingPts(){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain,iY;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator;
	CLog::Info("--- Y_train     Y_emulator    Sigma_emulator ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		for(iY=0;iY<NObservables;iY++){
			GetY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
			snprintf(pchars,CLog::CHARLENGTH,
			"Y[%u]=%10.3e =? %10.3e  +/- %12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
			CLog::Info(pchars);
		}
	}
}

void CSmoothMaster::TestAtTrainingPts(unsigned int iY){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain;
	double Y,SigmaY_emulator;
	CLog::Info("--- Y_train     Y_emulator    Sigma_emulator ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		GetY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e  +/- %12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
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
		GetY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::TestVsFullModelAlt(){
	char pchars[CLog::CHARLENGTH];
	unsigned int iY,ipar,nfit=0,ntest=0;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator,realY;
	vector<double> testtheta;
	FILE *fptr,*fptr_out;
	string filename;
	fitpercentage=0.0;
	for(iY=0;iY<NObservables;iY++){
		nfit=ntest=0;
		filename="smooth_data/fullmodel_testdata/"+observableinfo->observable_name[iY]+".txt";
		printf("filename=%s\n",filename.c_str());
		fptr=fopen(filename.c_str(),"r");
		filename="smooth_data/fullmodel_testdata/YvsY_"+observableinfo->observable_name[iY]+".txt";
		fptr_out=fopen(filename.c_str(),"w");
		
		testtheta.resize(NPars);
		do{
			for(ipar=0;ipar<NPars;ipar++){
				fscanf(fptr,"%lf",&testtheta[ipar]);
			}
			fscanf(fptr,"%lf",&realY);
			if(!feof(fptr)){
				ntest+=1;
				GetY(iY,testtheta,Y,SigmaY_emulator);
				snprintf(pchars,CLog::CHARLENGTH,
				"Y[%u]=%10.3e =? %10.3e,    SigmaY_emulator=%12.5e\n",
				iY,Y,realY,SigmaY_emulator);
				fprintf(fptr_out,"%12.5e  %12.5e %12.5e\n",realY,Y,SigmaY_emulator);
				if(fabs(Y-realY)<SigmaY_emulator)
					nfit+=1;
			}	
		}while(!feof(fptr));
		fclose(fptr);		
		fclose(fptr_out);
		CLog::Info(observableinfo->observable_name[iY]+": "+to_string(nfit)+" out of "+to_string(ntest)+" points within 1 sigma\n");
		fitpercentage+=100.0*double(nfit)/double(ntest);
	}
	fitpercentage=fitpercentage/double(NObservables);
	CLog::Info("percentage within 1 sigma = "+to_string(fitpercentage)+"\n");
}

void CSmoothMaster::TestVsFullModel(){
	string TestListStr = parmap->getS("SmoothEmulator_TestPts","1");
	
	vector<unsigned int> TestList;
	stringstream ss(TestListStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos) {

			unsigned int start = stoi(token.substr(0, pos));
			unsigned int end = stoi(token.substr(pos+1));

			for (unsigned int i = start; i <= end; i++)
			TestList.push_back(i);
		}
		else {
			TestList.push_back(stoi(token));
		}
	}
	
	unsigned int ntestpts=TestList.size();
	char obsnamechars[200],modparnamechars[200];
	string obsname,modparname;
	unsigned int iY,iread,itest,ipar,nfit;
	unsigned int NObservables=observableinfo->NObservables;
	double Y,SigmaY_emulator,realY,realSigmaY;
	double Xread,SigmaXRead;
	CModelParameters testpars;
	FILE *fptr,*fptr_out;
	string filename;
	for(iY=0;iY<NObservables;iY++){
		filename="smooth_data/fullmodel_testdata/YvsY_"+observableinfo->observable_name[iY]+".txt";
		fptr_out=fopen(filename.c_str(),"w");
		nfit=0;
		
		//CLog::Info("Writing test_vs_full_model results to "+filename+"\n");
		
		for(itest=0;itest<ntestpts;itest++){
			filename="smooth_data/modelruns/run"+to_string(TestList[itest])+"/mod_parameters.txt";
			fptr=fopen(filename.c_str(),"r");
			for(iread=0;iread<NPars;iread++){
				fscanf(fptr,"%s %lf %lf",modparnamechars,&Xread,&SigmaXRead);
				modparname=string(modparnamechars);
				ipar=priorinfo->GetIPosition(modparname);
				testpars.X[ipar]=Xread;
			}
			fclose(fptr);
			testpars.TranslateX_to_Theta();
			GetY(iY,testpars.Theta,Y,SigmaY_emulator);
			
			filename="smooth_data/modelruns/run"+to_string(TestList[itest])+"/obs.txt";
			fptr=fopen(filename.c_str(),"r");
			iread=-1;
			do{
				iread+=1;
				fscanf(fptr,"%s %lf %lf",obsnamechars,&realY,&realSigmaY);
				obsname=string(obsnamechars);
				
			}while(iread<observableinfo->NObservables && obsname!=observableinfo->observable_name[iY]);
			fclose(fptr);
			if(fabs(Y-realY)<SigmaY_emulator)
				nfit+=1;
			if(obsname!=observableinfo->observable_name[iY])
				CLog::Fatal("cannot find obsname amongst observales, obsname="+obsname+"\n");
			
			fprintf(fptr_out,"%lf %lf %lf\n",realY,Y,SigmaY_emulator);
			
		}
		fclose(fptr_out);
		CLog::Info(observableinfo->observable_name[iY]+": "+to_string(nfit)+" out of "+to_string(ntestpts)+" points within 1 sigma\n");
	}
}

vector<double> CSmoothMaster::GetThetaFromX(vector<double> X){
	double sigmax,xbar;
	unsigned int ipar;
	vector<double> Theta;
	Theta.resize(X.size());

	for(ipar=0;ipar<CModelParameters::NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			Theta[ipar]=-1+2*((X[ipar]-priorinfo->xmin[ipar])/(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]));
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			Theta[ipar]=(X[ipar]-xbar)/(sigmax*CModelParameters::GSCALE);
		}
		else{
			CLog::Fatal("Cannot translate X to Theta because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}
	return Theta;
		
}

vector<double> CSmoothMaster::GetXFromTheta(vector<double> Theta){
	double sigmax,xbar;
	unsigned int ipar;
	vector<double> X;
	X.resize(Theta.size());

	for(ipar=0;ipar<CModelParameters::NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			X[ipar]=priorinfo->xmin[ipar]+0.5*(1.0+Theta[ipar])*(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]);
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			X[ipar]=xbar+CModelParameters::GSCALE*sigmax*Theta[ipar];
		}
		else{
			CLog::Fatal("Cannot translate Theta to X because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}
	return X;
	
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



