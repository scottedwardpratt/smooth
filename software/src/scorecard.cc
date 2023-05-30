#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/scorecard.h"
using namespace std;


void CScoreCard::CalcScore(CSmoothEmulator *emulator,vector<vector<double>> &ThetaTest,double YExpSet,double SigmaYExpSet){
	YExp=YExpSet;
	SigmaYExp=SigmaYExpSet;
	int itest,isample,NTest;
	double yi,Pi,Pibar,Pi2bar;
	score=0.0;
	NTest=ThetaTest.size();
	if(ThetaTest.size()!=NTest)
		NTest=ThetaTest.size();
	
	for(itest=0;itest<NTest;itest++){
		Pibar=Pi2bar=0.0;
		for(isample=0;isample<emulator->NASample;isample++){
			yi=emulator->smooth->CalcY(emulator->ASample[isample],emulator->LAMBDA,ThetaTest[itest]);
			Pi=exp(-(yi-YExp)*(yi-YExp)/(2.0*SigmaYExp*SigmaYExp));
			Pibar+=Pi;
			Pi2bar+=Pi*Pi;
		}
		Pibar=Pibar/emulator->NASample;
		Pi2bar=Pi2bar/emulator->NASample;
		score+=sqrt(Pi2bar-Pibar*Pibar);
	}
}