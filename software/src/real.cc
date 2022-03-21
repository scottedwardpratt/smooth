#include "real.h"
#include "smooth.h"
#include "randy.h"
#include "emulator.h"
#include "constants.h"
#include "gslmatrix.h"

using namespace std;

CReal::CReal(){
//	cout << "CReal object created" << endl;
}

double CReal::CalcY(vector<double> &theta){
	cout << "object created2" << endl;
	return 0.0;
}

double CReal::CalcYReal(vector<double> &theta){
	cout << "object created2" << endl;
	return 0.0;
}

void CReal::RandomizeRealA(double SigmaReal){
	cout << "object created2" << endl;
//	return 0.0;
}

void CReal::CalcYTrainFromRealA(vector<double> YTrain,int NTrainingPts, vector<vector<double>> ThetaTrain){
	cout << "object created2" << endl;
//	return 0.0;
}

CReal_Taylor::CReal_Taylor(unsigned int NPars_Set,CRandy *randyset){	
//	cout << NPars << endl;
	NPars=NPars_Set;	
//	cout << "howdy, NPars=" << NPars << endl;
	randy=randyset;
	double crap=randy->ran_gauss();
//	cout << crap << endl;
//	cout << "object created2.5" << endl;
//	cout << NPars << endl;
	smooth = new CSmooth(NPars);
//	cout << "object created3" << endl;
//	RealA.resize(smooth->NCoefficients);
	LAMBDA=3;
/*	int NASample=10;
	RealA.resize(NASample);
	for(unsigned int isample=0;isample<NASample;isample++){
		ASample[isample].resize(smooth->NCoefficients);
		for(unsigned int ic=0;ic<ASample.size();ic++){
			A[ic]=0.0;
		}
	}
	*/
}


double CReal_Taylor::CalcYReal(vector<double> &theta){
//	cout << RealA.size() << endl;
	//size is 56
	
	return smooth->CalcY(RealA,LAMBDA,theta);
}

void CReal_Taylor::RandomizeRealA(double SigmaReal){
	if(RealA.size()!=smooth->NCoefficients){
		RealA.resize(smooth->NCoefficients);
	}
	for(unsigned int ic=0;ic<RealA.size();ic++){
		RealA[ic]=SigmaReal*randy->ran_gauss();
	}	
}

void CReal_Taylor::CalcYTrainFromRealA(vector<double> YTrain, int NTrainingPts, vector<vector<double>> ThetaTrain){
//	cout << "NTrainingPts" << NTrainingPts << endl;	
	//NtrainingPts is 4
	unsigned int iTrain;
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		YTrain[iTrain]=CalcYReal(ThetaTrain[iTrain]);
	}
}
