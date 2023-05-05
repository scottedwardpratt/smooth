#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/scorecard.h"
using namespace std;


void CScoreCard::CalcScore(CSmoothEmulator *emulator,CSmooth *smooth,vector<double> &ThetaTestSet,double YExpSet,double sigmaYExpSet){
	ThetaTest=ThetaTestSet;
	YExp=YExpSet;
	sigmaYExp=sigmaYExpSet;
	int itest,isample;
	double yi,Pi,Pibar,Pi2bar;
	score=0.0;
	if(ThetaTestSet.size()!=NTest)
		NTest=ThetaTestSet.size();
	
	for(itest=0;itest<NTest;itest++){
		Pibar=Pi2bar=0.0;
		for(isample=0;isample<emulator->NASample;isample++){
			yi=smooth->CalcY(emulator->ASample[isample],emulator->LAMBDA,ThetaTest);
			Pi=exp(-(yi-YExp)*(yi-YExp)/(2.0*sigmaYExp*sigmaYExp));
			Pibar+=Pi;
			Pi2bar+=Pi;
		}
		Pibar=Pibar/emulator->NASample;
		Pi2bar=Pi2bar/emulator->NASample;
		score+=sqrt(Pi2bar-Pibar*Pibar);
	}
}