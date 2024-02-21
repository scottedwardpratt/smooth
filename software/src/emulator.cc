#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;

unsigned int CSmoothEmulator::NPars=0;
CSmooth *CSmoothEmulator::smooth=NULL;
CSmoothMaster *CSmoothEmulator::smoothmaster=NULL;
CparameterMap *CSmoothEmulator::parmap=NULL;
Crandy *CSmoothEmulator::randy=NULL;
unsigned int CSmoothEmulator::NTrainingPts=0;

CSmoothEmulator::CSmoothEmulator(string observable_name_set){
	observable_name=observable_name_set;
	NTrainingPts=smoothmaster->traininginfo->NTrainingPts;

	LAMBDA=parmap->getD("SmoothEmulator_LAMBDA",3.0);
	NMC=parmap->getI("SmoothEmulator_MCMC_NMC",10000);
	NASample=parmap->getI("SmoothEmulator_MCMC_NASample",8);
	MCStepSize=parmap->getD("SmoothEmulator_MCStepSize",0.01);
	MCSigmaAStepSize=parmap->getD("SmoothEmulator_MCSigmaAStepSize",0.01);
	TuneChooseMCMC=parmap->getB("SmoothEmulator_TuneChooseMCMC",false);
	TuneChooseMCMCPerfect=parmap->getB("SmoothEmulator_TuneChooseMCMCPerfect",false);
	TuneChooseExact=parmap->getB("SmoothEmulator_TuneExact",true);
	UseSigmaY=parmap->getB("SmoothEmulator_MCMCUseSigmaY",false);
	ConstrainA0=parmap->getB("SmoothEmulator_ConstrainA0",false);
	CutOffA=parmap->getB("SmoothEmulator_MCMC_CutoffA",false);
	iY=smoothmaster->observableinfo->GetIPosition(observable_name);
	SigmaA0=smoothmaster->observableinfo->SigmaA0[iY];
	SigmaA=SigmaA0;
	SigmaAMin=parmap->getD("SmoothEmulator_SigmaAMin",0.1*SigmaA0);
	Init();
}

void CSmoothEmulator::SetThetaTrain(){
	ThetaTrain.resize(NTrainingPts);
	for(unsigned int itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(unsigned int ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=smoothmaster->traininginfo->modelpars[itrain]->Theta[ipar];
		}
	}
}

void CSmoothEmulator::Init(){
	SigmaA=SigmaA0;
	NSigmaA=0;
	SigmaAbar=0.0;
	FirstTune=true;
	MCStepSize=MCStepSize/double(NPars*NPars);
	MCSigmaAStepSize=MCSigmaAStepSize/double(NPars*NPars);

	ASample.resize(NASample);
	for(unsigned int isample=0;isample<NASample;isample++){
		ASample[isample].resize(smooth->NCoefficients);
		//SetA_RanGauss(SigmaA,ASample[isample]);
		SetA_Zero(ASample[isample]);
	}

	A.resize(smooth->NCoefficients);
	SetA_Zero(A);
	ATrial.resize(smooth->NCoefficients);
	SetA_Zero(ATrial);
	Ttilde.resize(NTrainingPts,NTrainingPts);
	TtildeInv.resize(NTrainingPts,NTrainingPts);
	T.resize(NTrainingPts);
	for(unsigned int it=0;it<NTrainingPts;it++){
		T[it].resize(smooth->NCoefficients);
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			T[it][ic]=0.0;
		}
	}
}

void CSmoothEmulator::Tune(){
	//FirstTune=true;
	CalcMForTraining();
	if(TuneChooseMCMC==true){
		if(UseSigmaY){
			if(FirstTune){
				TuneMCMC();
				FirstTune=false;
			}
			TuneMCMC_withSigma();
		}
		else{
			if(FirstTune){
				TuneMCMC();
				FirstTune=false;
			}
			TuneMCMC();
		}
	}
	else if(TuneChooseMCMCPerfect==true){
		TunePerfectMCMC();
	}
	else if(TuneChooseExact){
		TuneExact();
		GetExactSigmaA();
		CalcExactLogP();
	}
	else{
		CLog::Fatal("In CSmoothEmulator::Tune(), no tuning method specified\n");
	}
}

void CSmoothEmulator::TuneMCMC(){
	vector<double> *Aswitch,*Aptr,*ATrialptr;
	double dlp,r,SigmaAswitch;
	unsigned int success=0,ic,imc;
	double BestLogP,stepsize;
	for(ic=0;ic<smooth->NCoefficients;ic++){
		ATrial[ic]=A[ic];
	}
	CalcAFromTraining(A);
	Aptr=&A;
	ATrialptr=&ATrial;
	SigmaATrial=SigmaA0;
	double logP,logPTrial;
	logP=GetLog_AProb(*Aptr,SigmaA);
	BestLogP=-1000000.0;
	for(imc=0;imc<NMC;imc++){

		SigmaATrial=fabs(SigmaA+SigmaA0*MCSigmaAStepSize*randy->ran_gauss());
		if(SigmaATrial<SigmaAMin){
			SigmaATrial=SigmaA;
		}
		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			stepsize=MCStepSize*pow(LAMBDA,smooth->rank[ic]);
			(*ATrialptr)[ic]=SigmaATrial*( ((*Aptr)[ic]/SigmaA)+stepsize*randy->ran_gauss() );
		}
		//
		CalcAFromTraining(*ATrialptr);
		//Aptr=ATrialptr;
		
		//
		logPTrial=GetLog_AProb(*ATrialptr,SigmaATrial);
		dlp=logPTrial-logP;

		if(dlp>0.0){
			Aswitch=Aptr;
			Aptr=ATrialptr;
			ATrialptr=Aswitch;
			SigmaAswitch=SigmaA;
			SigmaA=SigmaATrial;
			SigmaATrial=SigmaAswitch;
			logP=logPTrial;
			success+=1;
		}
		else{
			if(dlp>-100){
				dlp=exp(dlp);
				r=randy->ran();
				if(dlp>r){
					Aswitch=Aptr;
					Aptr=ATrialptr;
					ATrialptr=Aswitch;
					SigmaAswitch=SigmaA;
					SigmaA=SigmaATrial;
					SigmaATrial=SigmaAswitch;
					logP=logPTrial;
					success+=1;
				}
			}
		}
		if(logP>BestLogP)
			BestLogP=logP;
	}

	if(Aptr!=&A){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=(*Aptr)[ic];
		}
	}

	unsigned int Ndof=smooth->NCoefficients-NTrainingPts;
	if(!FirstTune)
		CLog::Info("success percentage="+to_string(double(success)*100.0/double(NMC))+",SigmaA="+to_string(SigmaA)+",  logP/Ndof="+to_string(logP/double(Ndof))+", BestLogP/Ndof="+to_string(BestLogP/double(Ndof))+"\n");
}

void CSmoothEmulator::TuneMCMC_withSigma(){
	vector<double> *Aswitch,*Aptr,*ATrialptr;
	double dlp,r,SigmaAswitch,Y,stepsize;
	unsigned int success=0,ic,imc,itrain;
	double BestLogP;
	vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
	vector<double> SigmaYTrain=smoothmaster->traininginfo->SigmaYTrain[iY];
	for(ic=0;ic<smooth->NCoefficients;ic++){
		ATrial[ic]=A[ic];
	}
	Aptr=&A;
	ATrialptr=&ATrial;
	SigmaATrial=SigmaA;
	double logP,logPTrial;
	logP=GetLog_AProb(*Aptr,SigmaA);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		Y=smooth->CalcY_FromMtot(*Aptr,T[itrain]);
		logP-=0.5*(YTrain[itrain]-Y)*(YTrain[itrain]-Y)/(SigmaYTrain[itrain]*SigmaYTrain[itrain]);
	}

	BestLogP=-1000000000.0;
	for(imc=0;imc<NMC;imc++){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			stepsize=SigmaA*MCStepSize*pow(LAMBDA,smooth->rank[ic]);
			(*ATrialptr)[ic]=(*Aptr)[ic]+stepsize*randy->ran_gauss();
		}
		SigmaATrial=fabs(SigmaA+SigmaA0*MCSigmaAStepSize*randy->ran_gauss());
		for(ic=0;ic<smooth->NCoefficients;ic++){
			(*ATrialptr)[ic]*=(SigmaATrial/SigmaA);
		}
		logPTrial=GetLog_AProb(*ATrialptr,SigmaATrial);


		for(itrain=0;itrain<NTrainingPts;itrain++){
			Y=smooth->CalcY_FromMtot(*ATrialptr,T[itrain]);
			logPTrial-=0.5*(YTrain[itrain]-Y)*(YTrain[itrain]-Y)/(SigmaYTrain[itrain]*SigmaYTrain[itrain]);
		}

		dlp=logPTrial-logP;
		if(dlp>0.0){
			Aswitch=Aptr;
			Aptr=ATrialptr;
			ATrialptr=Aswitch;
			SigmaAswitch=SigmaA;
			SigmaA=SigmaATrial;
			SigmaATrial=SigmaAswitch;
			logP=logPTrial;
			success+=1;
		}
		else{
			if(dlp>-100){
				dlp=exp(dlp);
				r=randy->ran();
				if(dlp>r){
					Aswitch=Aptr;
					Aptr=ATrialptr;
					ATrialptr=Aswitch;
					SigmaAswitch=SigmaA;
					SigmaA=SigmaATrial;
					SigmaATrial=SigmaAswitch;
					logP=logPTrial;
					success+=1;
				}
			}
		}
		if(logP>BestLogP)
			BestLogP=logP;
	}

	if(Aptr!=&A){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=(*Aptr)[ic];
		}
	}

	unsigned int Ndof=smooth->NCoefficients-NTrainingPts;
	if(!FirstTune)
		CLog::Info("success percentage="+to_string(double(success)*100.0/double(NMC))+", SigmaA="+to_string(SigmaA)+", logP/Ndof="+to_string(logP/double(Ndof))+",BestLogP/Ndof="+to_string(BestLogP/double(Ndof))+"\n");
}

void CSmoothEmulator::TunePerfectMCMC(){
	unsigned int ic,ic0,ntry=0,ntrymax=100000;
	bool success=false;
	double weight,warg;//sigmafact=1.0;
	CalcMForTraining();
	if(ConstrainA0){
		ic0=0;
	}
	else
		ic0=1;

	while(ntry<ntrymax && success==false){
		SigmaA=SigmaA0;

		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			//SigmaA=0.5*SigmaA0*tan((PI/2.0)*(1.0-2.0*randy->ran()));
			ATrial[ic]=SigmaA*randy->ran_gauss();
		}
		warg=0.0;
		CalcAFromTraining(ATrial);
		for(ic=ic0;ic<NTrainingPts;ic++){
			warg-=0.5*ATrial[ic]*ATrial[ic]/(SigmaA*SigmaA);
		}
		//warg-=(NTrainingPts-1)*log(SigmaA/(0.5*SigmaA0));
		if(warg>0.0){
			CLog::Fatal("Disaster, warg="+to_string(warg)+"\n");
		}
		if(warg>-100){
			weight=exp(warg);
			if(weight>randy->ran()){
				success=true;
			}
		}
		ntry+=1;
	}
	if(ntry>=ntrymax){
		CLog::Fatal("TunePerfectMCMC Failed, SigmaA0="+to_string(SigmaA0)+"\n");
	}
	else{
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=ATrial[ic];
		}
	}
}

void CSmoothEmulator::CalcMForTraining(){
	unsigned int itrain,ic;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			T[itrain][ic]=smooth->GetT(ic,LAMBDA,ThetaTrain[itrain]);
		}
		for(ic=0;ic<NTrainingPts;ic++){
			Ttilde(itrain,ic)=T[itrain][ic];
		}
	}
	TtildeInv=Ttilde.inverse();
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ic=0;ic<NTrainingPts;ic++){
			if(TtildeInv(itrain,ic)!=TtildeInv(itrain,ic)){
				CLog::Fatal("TtildeInv != TtildeInv\n");
			}
		}
	}
}

// This adjust first NTrainingPts coefficients to reproduce Y training values

void CSmoothEmulator::CalcAFromTraining(vector<double> &AA){
	vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
	unsigned int itrain;
	//vector<double> Ashort;
	//Ashort.resize(smooth->NCoefficients);
	Eigen::VectorXd Ashort,YTarget;
	AA.resize(NTrainingPts);
	YTarget.resize(NTrainingPts);
	if(ThetaTrain.size()!=NTrainingPts){
		CLog::Info("CSmoothEmulator:: array size mismatch!!\n");
		CLog::Fatal("ThetaTrain.size="+to_string(ThetaTrain.size())+", NTrainingPts="+to_string(NTrainingPts)+"\n");
	}
	YTarget.resize(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		//YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,LAMBDA,ThetaTrain[itrain],NTrainingPts);
		YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder_FromMtot(AA,NTrainingPts,T[itrain]);
	}
	Ashort=TtildeInv*YTarget;
	for(itrain=0;itrain<NTrainingPts;itrain++)
		AA[itrain]=Ashort(itrain);
}

void CSmoothEmulator::OldCalcAFromTraining(vector<double> &AA){
	vector<double> YTrain=smoothmaster->traininginfo->YTrain[iY];
	unsigned int itrain,ic;
	//vector<double> Ashort;
	//Ashort.resize(smooth->NCoefficients);
	Eigen::VectorXd Ashort,YTarget;
	AA.resize(NTrainingPts);
	YTarget.resize(NTrainingPts);
	if(ThetaTrain.size()!=NTrainingPts){
		CLog::Info("CSmoothEmulator:: array size mismatch!!\n");
		CLog::Fatal("ThetaTrain.size="+to_string(ThetaTrain.size())+", NTrainingPts="+to_string(NTrainingPts)+"\n");
	}
	YTarget.resize(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,LAMBDA,ThetaTrain[itrain],NTrainingPts);
		for(ic=0;ic<NTrainingPts;ic++){
			Ttilde(itrain,ic)=smooth->GetT(ic,LAMBDA,ThetaTrain[itrain]);
		}
	}
	//Ashort=Ttilde.colPivHouseholderQr().solve(YTarget);
	Ashort=Ttilde.partialPivLu().solve(YTarget);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		AA[itrain]=Ashort(itrain);
}

double CSmoothEmulator::GetLog_AProb(vector<double> &AA,double ASigmaA){
	double answer=0.0;
	unsigned int ic0=1;
	if(ConstrainA0)
		ic0=0;
	//answer=0.000*smooth->NCoefficients*log(ASigmaA/SigmaA0);
	// Don't prefer small A[0], so start at ic=1
	for(unsigned int ic=ic0;ic<smooth->NCoefficients;ic++){
		answer-=0.5*AA[ic]*AA[ic]/(ASigmaA*ASigmaA);
	}
	if(!UseSigmaY){
		if(ConstrainA0)
			answer-=NTrainingPts*log(ASigmaA);
		else
			answer-=(NTrainingPts-1)*log(ASigmaA);
	}
	// next line keeps A from drifting out to infinity
	if(CutOffA)
			answer-=log(1.0+0.25*(ASigmaA*ASigmaA)/(SigmaA0*SigmaA0));
	return answer;
}

void CSmoothEmulator::SetA_Zero(vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=0.0;
}

void CSmoothEmulator::SetA_RanGauss(double SigmaA,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaA*randy->ran_gauss();
}

void CSmoothEmulator::SetA_Constant(double SigmaA,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaA;
}

void CSmoothEmulator::SetA_RanSech(double SigmaA,vector<double> &A){
	double r=1.0-2.0*randy->ran();
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaA*2.0*atanh(tan(0.5*PI*(r-0.5)));
}

void CSmoothEmulator::GenerateASamples(){
	unsigned int isample;
	FirstTune=true;
	for(isample=0;isample<NASample;isample++){
		Tune();
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			ASample[isample][ic]=A[ic];
		}
	}
	SigmaAbar+=SigmaA;
	NSigmaA+=1;
}

