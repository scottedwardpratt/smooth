#include "emulator.h"
#include "smooth.h"
using namespace std;

double CSmoothEmulator::CalcRealYFromRealA(vector<double> &theta){
	return smooth->CalcY(RealA,LAMBDA,theta);
}

void CSmoothEmulator::RandomizeRealA(){
if(RealA.size()!=smooth->NCoefficients)
		RealA.resize(smooth->NCoefficients);
	SetA_RanGauss(SigmaY0,RealA);
}

void CSmoothEmulator::CalcYTrainFromRealA(){
	unsigned int iTrain;
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		YTrain[iTrain]=CalcRealYFromRealA(ThetaTrain[iTrain]);
	}
}