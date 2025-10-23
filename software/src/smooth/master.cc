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
   string filename;
   
   filename="smooth_data/Info/observable_info.txt";
   
   observableinfo=new CObservableInfo(filename);
   NObs=observableinfo->NObservables;
   
	FullModelRunDirName=parmap->getS("Smooth_FullModelRunDirName","smooth_data/FullModelRuns");
   FullModelTestingRunDirName=parmap->getS("Smooth_FullModelTestingRunDirName","smooth_data/FullModelTestingRuns");
   
	SurmiseTrainingPointsFileName=parmap->getS("SmoothEmulator_SurmiseTrainingPointsFilename","SurmiseTrainingPoints.txt");
	SurmiseTrainingObsFileName=parmap->getS("SmoothEmulator_SurmiseTrainingObsFilename","SurmiseTrainingObs.txt");
   
   SurmiseTestingPointsFileName=parmap->getS("SmoothEmulator_SurmiseTestingPointsFilename","SurmiseTestingingPoints.txt");
   SurmiseTestingObsFileName=parmap->getS("SmoothEmulator_SurmiseTestingObsFilename","SurmiseTestingingObs.txt");
	
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
	FILE *fptr=fopen("sigmalambda.txt","w");
	double sigmaAbar=0.0,Lambdabar=0.0;
	for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
		emulator[iY]->Tune();
		fprintf(fptr,"%10.3f %10.5f\n",emulator[iY]->SigmaA,emulator[iY]->LAMBDA);
		//printf("%10.3f %10.5f\n",emulator[iY]->SigmaA,emulator[iY]->LAMBDA);
		sigmaAbar+=emulator[iY]->SigmaA;
		Lambdabar+=emulator[iY]->LAMBDA;
	}
	sigmaAbar=sigmaAbar/double(observableinfo->NObservables);
	Lambdabar=Lambdabar/double(observableinfo->NObservables);
	fclose(fptr);
	//printf("<sigmaA>=%g, <Lambda>=%g\n",sigmaAbar,Lambdabar);
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

void CSmoothMaster::TestAtTrainingPts(){
	unsigned int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		TestAtTrainingPts(iY);
	}

}

void CSmoothMaster::TestAtTrainingPts(unsigned int iY){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain;
	double Y,SigmaY_emulator;
	CLog::Info("----------- "+observableinfo->observable_name[iY]+" -------------\n");
	CLog::Info("--- Y_train     Y_emulator    Sigma_emulator ---- Lambda="
		+to_string(emulator[iY]->LAMBDA)+", SigmaA="+to_string(emulator[iY]->SigmaA)+"\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+": ");
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

void CSmoothMaster::TestVsFullModel(){
   char obsnameread[200],cdummy[200];
	unsigned int iY,ntest=0,itest;
	unsigned int NObs=observableinfo->NObservables;
	vector<double> Y,SigmaY_emulator;
	vector<double> testtheta;
	FILE *fptr,*fptr_out;
	string filename;
	fitpercentage=0.0;
   vector<int> nfit;
   vector<double> averagesigma2;
   string(command);
   double realY;
	
	int ifit;
	int nfitpercent[100]={0};
   ReadTestingInfo();
   ntest=CSmoothEmulator::NTestingPts;
   nfit.resize(NObs,0);
   averagesigma2.resize(NObs,0.0);
   Y.resize(NObs);
   SigmaY_emulator.resize(NObs);
   
   command="mkdir -p smooth_data/fullmodel_testdata";
   filename=FullModelTestingRunDirName+"/YvsY.txt";
   fptr_out=fopen(filename.c_str(),"w");
   for(itest=0;itest<ntest;itest++){
      GetAllY(testinginfo->modelpars[itest],Y,SigmaY_emulator); // emulated values
      
      
      filename=FullModelTestingRunDirName+"/run"+to_string(itest)+"/obs.txt";
      fptr=fopen(filename.c_str(),"r");
      do{
         fscanf(fptr,"%s",obsnameread);
         while(obsnameread[0]=='#'){
            fgets(cdummy,200,fptr);
            fscanf(fptr,"%s",obsnameread);
         }
         if(!feof(fptr)){
            fscanf(fptr,"%lf",&realY);
            if(!feof(fptr)){
               iY=observableinfo->GetIPosition(obsnameread);
               averagesigma2[iY]+=(realY-Y[iY])*(realY-Y[iY]);
               fprintf(fptr_out,"%10.3e %10.3e %10.3e\n",realY,Y[iY],SigmaY_emulator[iY]);
               if(fabs(realY-Y[iY])<SigmaY_emulator[iY])
                  nfit[iY]+=1;
            }
            fgets(cdummy,200,fptr);
            fscanf(fptr,"%s",obsnameread);
         }
         
         
         
      }while(!feof(fptr));
      fclose(fptr);
      
   }
   fclose(fptr_out);
   
   
   for(iY=0;iY<NObs;iY++){
      CLog::Info("<(Y-Yreal)^2>/SigmaY^2_emulator="+to_string(averagesigma2[iY]/(SigmaY_emulator[iY]*SigmaY_emulator[iY]))+"\n");
      fitpercentage=100.0*double(nfit[iY])/double(ntest);
      CLog::Info("percent < 1 sigma = "+to_string(fitpercentage)+"\n");
      ifit=floorl(fitpercentage);
      nfitpercent[ifit]+=1;
   }
   
   filename="fitpercentage.txt";
   fptr=fopen(filename.c_str(),"w");
   for(ifit=0;ifit<100;ifit++){
      fprintf(fptr,"%3d %g\n",ifit,double(nfitpercent[ifit])/double(ntest));
   }
   
   fclose(fptr);

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



