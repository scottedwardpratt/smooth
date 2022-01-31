#include "real.h"
#include "emulator.h"
#include "smooth.h"

using namespace std;

CReal::CReal(){
	cout << "object created" << endl;
//	CSmooth *smooth;
}

CReal_Taylor::CReal_Taylor(unsigned int NPars_Set){
	NPars=NPars_Set;
	CSmooth *smooth = new CSmooth(NPars);
	
}
double CSmoothEmulator::CalcRealYFromRealA(vector<double> &theta){
	return smooth->CalcY(RealA,LAMBDA,theta);
}

void CSmoothEmulator::CalcYTrainFromRealA(){
	unsigned int iTrain;
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		YTrain[iTrain]=CalcRealYFromRealA(ThetaTrain[iTrain]);
	}
}

void CSmoothEmulator::RandomizeRealA(){
	if(RealA.size()!=smooth->NCoefficients)
			RealA.resize(smooth->NCoefficients);
		SetA_RanGauss(SigmaY0,RealA);
}

/*
double CSmoothEmulator::CalcRealYFromRealA(vector<double> &theta){
	return smooth->CalcY(RealA,LAMBDA,theta);
}


void CSmoothEmulator::CalcYTrainFromRealA(){
	unsigned int iTrain;
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		YTrain[iTrain]=CalcRealYFromRealA(ThetaTrain[iTrain]);
	}
}

void CSmoothEmulator::RandomizeRealA(){
if(RealA.size()!=smooth->NCoefficients)
		RealA.resize(smooth->NCoefficients);
	SetA_RanGauss(SigmaY0,RealA);
}
*/ 