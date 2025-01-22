#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::CalcSigmaA(){
	CalcB();
	int a,b;
	double sigmaA2=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		for(b=0;b<int(NTrainingPts);b++){
			sigmaA2+=smoothmaster->traininginfo->YTrain[iY][a]
				*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	sigmaA2=sigmaA2/double(NTrainingPts);
	SigmaA=sqrt(fabs(sigmaA2));

}

void CSmoothEmulator::CalcLogP(){
	double exparg=0.0;
	for(int a=0;a<int(NTrainingPts);a++){
		for(int b=0;b<int(NTrainingPts);b++){
			exparg-=0.5*smoothmaster->traininginfo->YTrain[iY][a]*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	exparg=exparg/(SigmaA*SigmaA);
	double detB=B.determinant();
	logP=-0.5*log(fabs(detB))-NTrainingPts*log(SigmaA)+exparg;
	if(isinf(logP)){
		printf("XXXXXXXXX detB=%g, Lambda=%g, SigmaA=%g, exparg=%g\n",detB,LAMBDA,SigmaA,exparg);
	}
}

void CSmoothEmulator::CalcSigmaALambda(){
	double LambdaMin=0.5*sqrt(double(NPars));
	if(LambdaMin<1.0)
		LambdaMin=1.0;
	double bestLambda,dLambda=1.0,bestlogP,oldbestlogP,oldbestLambda;
	int nfail=0;
	LAMBDA=LambdaMin;  // minimum LAMBDA
	CalcB();
	CalcSigmaA();
	CalcLogP();
	bestLambda=LAMBDA;
	bestlogP=logP;
	oldbestlogP=logP;
	oldbestLambda=LAMBDA;
	if(logP!=logP || isinf(logP)){
		nfail=100;
		CLog::Info("Note: LAMBDA set to "+to_string(LAMBDA)+" -- optimum value might be higher\n but determinant(B) cannot be calculated due to numerical accuracy problems.\n");
		logP=-200.0;
	}
	LAMBDA=bestLambda+dLambda;
	while(nfail<7){
		CalcB();
		CalcSigmaA();
		CalcLogP();
		if(logP==logP && !isinf(logP) && logP>bestlogP){
			oldbestlogP=bestlogP;
			oldbestLambda=bestLambda;
			bestlogP=logP;
			bestLambda=LAMBDA;
			LAMBDA=bestLambda+dLambda;
		}
		else{
			bestLambda=oldbestLambda;
			bestlogP=oldbestlogP;
			nfail+=1;
			dLambda=dLambda/2.0;
			LAMBDA=bestLambda+dLambda;
			if(logP!=logP || isinf(logP)){
				CLog::Info("Note: LAMBDA set to "+to_string(LAMBDA)+" -- optimum value might be higher\n but determinant(B) cannot be calculated due to numerical accuracy problems.\n");
			}
		}
	}
	LAMBDA=bestLambda;
	CalcB();
	CalcSigmaA();
	CalcLogP();
}

void CSmoothEmulator::CalcLambdaVariance(){
	double L,Lbar,dL=0.1,L2bar,norm,w;
	norm=Lbar=L2bar=0.0;
	for(L=1.0;L<8;L+=dL){
		LAMBDA=L;
		CalcB();
		CalcSigmaA();
		CalcLogP();
		//w=exp(logP);
		w=sqrt(Binv.determinant());
		double factor=pow(SigmaA/100.0,NTrainingPts);
		w=w/factor;
		printf("L=%g, w=%g, |Binv|=%g, SigmaA=%g\n",L,w,Binv.determinant(),SigmaA);
		norm+=w;
		Lbar+=w*L;
		L2bar+=w*L*L;
	}
	Lbar=Lbar/norm;
	L2bar=L2bar/norm;
	LambdaVariance=L2bar-Lbar*Lbar;
	LAMBDA=Lbar;
}
