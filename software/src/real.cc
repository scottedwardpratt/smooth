#include "real.h"
#include "smooth.h"

using namespace std;

CReal::CReal(){
	cout << "object created" << endl;
//	CSmooth *smooth;
}

CReal_Taylor::CReal_Taylor(unsigned int NPars_Set,CRandy *randyset){
	NPars=NPars_Set;
	randy=randyset;
	smooth = new CSmooth(NPars);
	RealA.resize(smooth->NCoefficients);
}

double CReal_Taylor::CalcY(vector<double> &theta){
	return smooth->CalcY(RealA,LAMBDA,theta);
}

void CReal_Taylor::RandomizeRealA(double SigmaReal){
	if(RealA.size()!=smooth->NCoefficients)
		RealA.resize(smooth->NCoefficients);
	for(unsigned int ic=0;ic<RealA.size();ic++)
		RealA[ic]=SigmaReal*randy->ran_gauss();
}
