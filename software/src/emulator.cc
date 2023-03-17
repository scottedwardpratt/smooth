#include "emulator.h"
#include "smooth.h"
using namespace std;

CSmoothEmulator::CSmoothEmulator(CparameterMap *parmap){
	NPars=parmap->getD("SmoothEmulator_NPars",0);
	parmap->set("Smooth_NPars",NPars);
	LAMBDA=parmap->getD("SmoothEmulator_LAMBDA",1.0);
	NMC=parmap->getI("SmoothEmulator_NMC",10000);
	NASample=parmap->getI("SmoothEmulator_NASample",10);
	MCStepSize=parmap->getD("SmoothEmulator_MCStepSize",0.5);
	MCSigmaYStepSize=parmap->getD("SmoothEmulator_MCSigmaYStepSize",0.1);
	SigmaY0=parmap->getD("SmoothEmulator_SigmaY",1.0);
	TuneAChooseMCMC=parmap->getB("SmoothEmulator_TuneAChooseMCMC",true);
	ConstrainA0=parmap->getB("SmoothEmulator_ConstrainA0",false);
	CutOffA=parmap->getB("SmoothEmulator_CutoffA",false);
	SigmaYMin=parmap->getD("SmoothEmulator_SigmaYMin",0.1*SigmaY0);

	smooth=new CSmooth(parmap);
	randy=new Crandy(-time(NULL));
	SigmaY=SigmaY0;
	NSigmaY=0;
	SigmaYbar=0.0;
	MCStepSize=MCStepSize/double(NPars*NPars);
	MCSigmaYStepSize=MCSigmaYStepSize/double(NPars*NPars);
//	cout << "sigmaY: " << SigmaY << endl;

	ASample.resize(NASample);
	simplex=new CSimplexSampler(parmap);
	
	for(unsigned int isample=0;isample<NASample;isample++){		
		ASample[isample].resize(smooth->NCoefficients);
		//SetA_RanGauss(SigmaY,ASample[isample]);
		SetA_Zero(ASample[isample]);		
	}
	
	A.resize(smooth->NCoefficients);
	//SetA_RanGauss(SigmaY,A);
	SetA_Zero(A);
	ATrial.resize(smooth->NCoefficients);
	SetA_Zero(ATrial);
}

void CSmoothEmulator::Init(CSmooth *smooth){
	randy=new Crandy(-time(NULL));
	SigmaY=SigmaY0;
	NSigmaY=0;
	SigmaYbar=0.0;
	MCStepSize=MCStepSize/double(NPars*NPars);
	MCSigmaYStepSize=MCSigmaYStepSize/double(NPars*NPars);

	ASample.resize(NASample);
//	simplex=new CSimplexSampler(parmap);
	for(unsigned int isample=0;isample<NASample;isample++){
		ASample[isample].resize(smooth->NCoefficients);
		//SetA_RanGauss(SigmaY,ASample[isample]);
		SetA_Zero(ASample[isample]);
	}
	A.resize(smooth->NCoefficients);
	//SetA_RanGauss(SigmaY,A);
	SetA_Zero(A);
	ATrial.resize(smooth->NCoefficients);
	SetA_Zero(ATrial);
	real=NULL;
}

void CSmoothEmulator::SetNTrainingPts(unsigned int NTrainingPts_set){
	NTrainingPts=NTrainingPts_set;
	M.resize(NTrainingPts,NTrainingPts);
	YTrain.resize(NTrainingPts);
	ThetaTrain.resize(NTrainingPts);
	for(unsigned int itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(unsigned int ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain].resize(NPars);
	}
}

void CSmoothEmulator::SetThetaSimplex(){
	unsigned int NTrain;
	simplex->SetThetaSimplex(ThetaTrain,NTrain);
	SetNTrainingPts(NTrain);
}

void CSmoothEmulator::TuneA(){
	if(TuneAChooseMCMC==true){
		TuneAMCMC();
	}
	else{
		TuneAPerfect();
	}
}

void CSmoothEmulator::TuneAMCMC(){
	vector<double> *Aswitch,*Aptr,*ATrialptr;
	double dlp,r,SigmaYswitch;
	unsigned int success=0,ic,imc;
	double BestLogP;
	CalcAFromTraining(A);
	for(ic=0;ic<smooth->NCoefficients;ic++){
		ATrial[ic]=A[ic];
	}
	Aptr=&A;
	ATrialptr=&ATrial;
	SigmaYTrial=SigmaY;
	double logP,logPTrial;
	logP=GetLog_AProb(*Aptr,SigmaY);
	BestLogP=logP;
	for(imc=0;imc<NMC;imc++){
		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			(*ATrialptr)[ic]=(*Aptr)[ic]+SigmaY*MCStepSize*randy->ran_gauss();
		}
		SigmaYTrial=fabs(SigmaY+SigmaY0*MCSigmaYStepSize*randy->ran_gauss());
		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			(*ATrialptr)[ic]*=(SigmaYTrial/SigmaY);
		}
		CalcAFromTraining(*ATrialptr);
		logPTrial=GetLog_AProb(*ATrialptr,SigmaYTrial);
		dlp=logPTrial-logP;	
			
		if(dlp>0.0){
			Aswitch=Aptr;
			Aptr=ATrialptr;
			ATrialptr=Aswitch;
			SigmaYswitch=SigmaY;
			SigmaY=SigmaYTrial;
			SigmaYTrial=SigmaYswitch;
			logP=logPTrial;
			success+=1;
//			cout << "dlp is:" << dlp << endl;
		}
		else{
			if(dlp>-100){
				dlp=exp(dlp);
				r=randy->ran();
//				cout << "dlp is:" << dlp << endl;
//				cout << "randy is: " << r <<endl;
				if(dlp>r){
					Aswitch=Aptr;
					Aptr=ATrialptr;
					ATrialptr=Aswitch;
					SigmaYswitch=SigmaY;
					SigmaY=SigmaYTrial;
					SigmaYTrial=SigmaYswitch;
					logP=logPTrial;
					success+=1;
				}
			}
		}
		if(logP>BestLogP)
			BestLogP=logP;
	}
//	cout << "dlp is now:" << dlp << endl;
	
	if(Aptr!=&A){
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=(*Aptr)[ic];
		}
	}
	
	printf("success percentage=%g, SigmaY=%g, logP=%g\n",double(success)*100.0/double(NMC),SigmaY,logP);
}

void CSmoothEmulator::TuneAPerfect(){
	unsigned int ic,ic0,Ns,ntry=0,ntrymax=1000000;
	bool success=false;
	double weight,warg;//sigmafact=1.0;
	if(ConstrainA0){
		ic0=0;
	}
	else
		ic0=1;
	Ns=NTrainingPts-1-ic0;

	while(ntry<ntrymax && success==false){
		if(Ns==0)
			SigmaY=SigmaY0;
		else
			SigmaY=SigmaYMin*pow(randy->ran(),-1.0/double(Ns));
		//SigmaY=0.5*SigmaY0*(1.0+2.0*randy->ran());

		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			//SigmaY=0.5*SigmaY0*tan((PI/2.0)*(1.0-2.0*randy->ran()));
			ATrial[ic]=SigmaY*randy->ran_gauss();
		}
		warg=0.0;
		CalcAFromTraining(ATrial);
		for(ic=ic0;ic<NTrainingPts;ic++){
			warg-=0.5*ATrial[ic]*ATrial[ic]/(SigmaY*SigmaY);
		}
		//warg-=(NTrainingPts-1)*log(SigmaY/(0.5*SigmaY0));
		if(warg>0.0){
			printf("Disaster, warg=%g\n",warg);
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
		printf("TuneAPerfect Failed, SigmaY0=%g\n",SigmaY0);
		exit(1);
	}
	else{
		for(ic=0;ic<smooth->NCoefficients;ic++){
			A[ic]=ATrial[ic];
		}
	}
}

// This adjust first NTrainingPts coefficients to reproduce Y training values
void CSmoothEmulator::CalcAFromTraining(vector<double> &AA){
	unsigned int itrain,ic;
	//vector<double> Ashort;
	//Ashort.resize(smooth->NCoefficients);
	Eigen::VectorXd Ashort,YTarget;
	AA.resize(NTrainingPts);
	YTarget.resize(NTrainingPts);
	if(ThetaTrain.size()!=NTrainingPts){
		printf("CSmoothEmulator:: array size mismatch!!\n");
		printf("ThetaTrain.size=%lu, NTrainingPts=%u\n",ThetaTrain.size(),NTrainingPts);
		exit(1);
	}
	YTarget.resize(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,LAMBDA,ThetaTrain[itrain],NTrainingPts);
		for(ic=0;ic<NTrainingPts;ic++){
			M(itrain,ic)=smooth->GetM(ic,LAMBDA,ThetaTrain[itrain]);
		}
	}
	Ashort=M.colPivHouseholderQr().solve(YTarget);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		AA[itrain]=Ashort(itrain);
}

double CSmoothEmulator::GetLog_AProb(vector<double> &AA,double ASigmaY){
	double answer=0.0;
	// Don't prefer small A[0], so start at ic=1
	for(unsigned int ic=1;ic<smooth->NCoefficients;ic++){
		//answer-=log(1.0+A[ic]*A[ic]/(SigmaY*SigmaY));
		answer-=0.5*AA[ic]*AA[ic]/(ASigmaY*ASigmaY);
	}
	answer-=(NTrainingPts-1)*log(ASigmaY);
	// next line keeps A from drifting out to infinity
	if(CutOffA)
		answer-=log(1.0+0.25*(ASigmaY*ASigmaY)/(SigmaY0*SigmaY0));
	return answer;
}

void CSmoothEmulator::SetA_Zero(vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=0.0;
}

void CSmoothEmulator::SetA_RanGauss(double SigmaY,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaY*randy->ran_gauss();
}

void CSmoothEmulator::SetA_Constant(double SigmaY,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaY;
}

void CSmoothEmulator::SetA_RanSech(double SigmaY,vector<double> &A){
	double r=1.0-2.0*randy->ran();
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=SigmaY*2.0*atanh(tan(0.5*PI*(r-0.5)));
}

void CSmoothEmulator::GenerateASamples(){
	unsigned int isample;
//	cout << "NASample is:" << NASample << endl;
	//NASample = 10
	for(isample=0;isample<NASample;isample++){
		TuneA();
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			ASample[isample][ic]=A[ic];
		}
	}
	SigmaYbar+=SigmaY;
	NSigmaY+=1;
}

void CSmoothEmulator::PrintA(vector<double> &Aprint){
	for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
		printf("%3u %g\n",ic,Aprint[ic]);
	}
}

void CSmoothEmulator::CalcYTrainFromThetaTrain(){
	unsigned int itrain;
	if(real==NULL){
		printf("Reality does not exist\n");
		exit(1);
	}
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTrain[itrain]=real->CalcY(ThetaTrain[itrain]);
	}
}