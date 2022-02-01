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
	smooth = new CSmooth(NPars);
	
}

double CReal_Taylor::CalcRealY(vector<double> &A, vector<double> &theta){
	double LAMBDA=3.0;
	return smooth->CalcY(A,LAMBDA,theta);
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