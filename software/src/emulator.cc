#include "emulator.h"
#include "smooth.h"
using namespace std;

CSmoothEmulator::CSmoothEmulator(CparameterMap *parmap){
	NPars=parmap->getD("Smooth_NPars",0);
	LAMBDA=parmap->getD("Smooth_LAMBDA",1.0);
	smooth=new CSmooth(NPars);
	randy=new CRandy(-time(NULL));
	MCStepSize=parmap->getD("Smooth_MCStepSize",0.5);
	MCStepSize=MCStepSize/double(NPars);
	MCAmagStepSize=parmap->getD("Smooth_MCAmagStepSize",0.1);
	MCAmagStepSize=MCAmagStepSize/double(NPars);
	Amag0=parmap->getD("Smooth_Amag",1.0);
	Amag=Amag0;
	NMC=parmap->getI("Smooth_NMC",10000);
	NASample=parmap->getI("Smooth_NASample",10);
	TrainRank=parmap->getI("Smooth_TrainRank",1);
	RTrain=parmap->getD("Smooth_RTrain",0.9);  
	ASample.resize(NASample);
	Lambda.resize(NPars);
	simplex=new CSimplexSampler(parmap);
	SetLambda_Constant(LAMBDA);
	for(unsigned int isample=0;isample<NASample;isample++){
		ASample[isample].resize(smooth->NCoefficients);
		//SetA_RanGauss(Amag,ASample[isample]);
		SetA_Zero(ASample[isample]);
	}
	A.resize(smooth->NCoefficients);
	//SetA_RanGauss(Amag,A);
	SetA_Zero(A);
	ATrial.resize(smooth->NCoefficients);
	SetA_Zero(ATrial);
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
	if(TrainRank==1)
		simplex->SetThetaRank1(ThetaTrain,NTrain);
	else if(TrainRank==2)
		simplex->SetThetaRank2(ThetaTrain,NTrain);
	else{
		printf("Smooth_ThetaRank in parmap must be 0 or 1\n");
		exit(1);
	}
	SetNTrainingPts(NTrain);
}

void CSmoothEmulator::TuneA(){
	vector<double> *Aswitch,*Aptr,*ATrialptr;
	double dlp,r,Amagswitch;
	unsigned int success=0,ic,imc;
	double BestLogP;
	CalcAFromTraining(A);
	for(ic=0;ic<smooth->NCoefficients;ic++){
		ATrial[ic]=A[ic];
	}
	Aptr=&A;
	ATrialptr=&ATrial;
	AmagTrial=Amag;
	double logP,logPTrial;
	logP=GetLog_AProb(*Aptr,Amag);
	BestLogP=logP;
	for(imc=0;imc<NMC;imc++){
		for(ic=NTrainingPts;ic<smooth->NCoefficients;ic++){
			(*ATrialptr)[ic]=(*Aptr)[ic]+MCStepSize*randy->ran_gauss();
		}
		//AmagTrial=sqrt(fabs(Amag*Amag+MCStepSize*randy->ran_gauss()));
		AmagTrial=fabs(Amag+Amag0*MCAmagStepSize*randy->ran_gauss());
		CalcAFromTraining(*ATrialptr);
		logPTrial=GetLog_AProb(*ATrialptr,AmagTrial);
		dlp=logPTrial-logP;
		if(dlp>0.0){
			Aswitch=Aptr;
			Aptr=ATrialptr;
			ATrialptr=Aswitch;
			Amagswitch=Amag;
			Amag=AmagTrial;
			AmagTrial=Amagswitch;
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
					Amagswitch=Amag;
					Amag=AmagTrial;
					AmagTrial=Amagswitch;
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
	
	printf("success percentage=%g, Amag=%g, logP=%g\n",double(success)*100.0/double(NMC),Amag,logP);
}

// This adjust first NTrainingPts coefficients to reproduce Y training values
void CSmoothEmulator::CalcAFromTraining(vector<double> &AA){
	unsigned int itrain,ic,ir;
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
		//YTarget[itrain]=YTrain[itrain]-smooth->CalcY_Remainder(AA,Lambda,ThetaTrain[itrain],NTrainingPts);
		YTarget(itrain)=YTrain[itrain]-smooth->CalcY_Remainder(AA,Lambda,ThetaTrain[itrain],NTrainingPts);
		for(ic=0;ic<NTrainingPts;ic++){
			//M[itrain][ic]=smooth->dupfactor[ic]/double(smooth->factorial[smooth->rank[ic]]);
			M(itrain,ic)=smooth->dupfactor[ic]/double(smooth->factorial[smooth->rank[ic]]);
			for(ir=0;ir<smooth->rank[ic];ir++){
				//M[itrain][ic]*=ThetaTrain[itrain][smooth->IPar[ic][ir]]/Lambda[smooth->IPar[ic][ir]];
				M(itrain,ic)*=ThetaTrain[itrain][smooth->IPar[ic][ir]]/Lambda[smooth->IPar[ic][ir]];
			}
		}
	}
	Ashort=M.colPivHouseholderQr().solve(YTarget);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		AA[itrain]=Ashort(itrain);
}

double CSmoothEmulator::GetLog_AProb(vector<double> &A,double Amag){
	double answer=0.0;
	// Don't prefer small A[0], so start at ic=1
	for(unsigned int ic=1;ic<smooth->NCoefficients;ic++){
		//answer-=log(1.0+A[ic]*A[ic]/(Amag*Amag));
		answer-=0.5*A[ic]*A[ic]/(Amag*Amag);
	}
	answer-=(smooth->NCoefficients-1)*log(Amag);
	answer-=Amag/Amag0;
	return answer;
}

void CSmoothEmulator::SetLambda_Constant(double LAMBDA_set){
	for(unsigned int ipar=0;ipar<NPars;ipar++)
		Lambda[ipar]=LAMBDA_set;
}

void CSmoothEmulator::SetA_Zero(vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=0.0;
}

void CSmoothEmulator::SetA_RanGauss(double Amag,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=Amag*randy->ran_gauss();
}

void CSmoothEmulator::SetA_Constant(double Amag,vector<double> &A){
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=Amag;
}

void CSmoothEmulator::SetA_RanSech(double Amag,vector<double> &A){
	double r=1.0-2.0*randy->ran();
	for(unsigned int ic=0;ic<A.size();ic++)
		A[ic]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
}

void CSmoothEmulator::GenerateASamples(){
	unsigned int isample;
	for(isample=0;isample<NASample;isample++){
		TuneA();
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			ASample[isample][ic]=A[ic];
		}
	}
}

void CSmoothEmulator::PrintA(vector<double> &Aprint){
	for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
		printf("%3u %g\n",ic,Aprint[ic]);
	}
}