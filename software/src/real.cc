#include "real.h"
#include "smooth.h"
#include "msu_commonutils/randy.h"
#include "emulator.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/gslmatrix.h"

using namespace std;

CReal::CReal(){
//	cout << "CReal object created" << endl;
}

double CReal::CalcY(vector<double> &theta){
	cout << "dummy function -- should not be hear" << endl;
	return 0.0;
}

CReal_Taylor::CReal_Taylor(unsigned int NPars_Set,Crandy *randyset){	
	NPars=NPars_Set;	
	randy=randyset;
	smooth = new CSmooth(NPars);
	LAMBDA=2;
}


double CReal_Taylor::CalcY(vector<double> &theta){
	return smooth->CalcY(A,LAMBDA,theta);
}

void CReal_Taylor::RandomizeA(double SigmaReal){
	if(A.size()!=smooth->NCoefficients){
		A.resize(smooth->NCoefficients);
	}
	for(unsigned int ic=0;ic<A.size();ic++){
		A[ic]=SigmaReal*randy->ran_gauss();
	}
}

void CReal::CalcYTrain(vector<double> &YTrain, int NTrainingPts, vector<vector<double>> ThetaTrain){
//	cout << "NTrainingPts" << NTrainingPts << endl;	
	//NtrainingPts is 4
	unsigned int itrain;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTrain[itrain]=CalcY(ThetaTrain[itrain]);
	}
}
