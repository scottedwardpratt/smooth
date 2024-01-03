#include "msu_smooth/smooth.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

double CSmooth::CalcY(vector<double> &A,double LAMBDA,vector<double> &theta){
	unsigned int ic,ir;
	double answer=0.0,term,rfactor;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=0.0;
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
//			cout << "theta is:" << theta[IPar[ic][ir]] << endl;
		}
		answer+=term;
	}
	answer*=rfactor;


	return answer;
}

double CSmooth::CalcY_Remainder(vector<double> &A,double LAMBDA,vector<double> &theta,unsigned int NTrainingPts){
	unsigned int ic,ir;
	double answer=0.0,term,rfactor;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=0.0;
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
		}
		answer+=term;
	}
	answer*=rfactor;
	return answer;
}

double CSmooth::CalcY_Remainder_FromMtot(vector<double> &A,unsigned int NTrainingPts,vector<double> &Mtot){
	double answer=0.0;
	unsigned int ic;
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		answer+=Mtot[ic]*A[ic];
	}
	return answer;
}

double CSmooth::CalcY_FromMtot(vector<double> &A,vector<double> &Mtot){
	double answer=0.0;
	unsigned int ic;
	for(ic=0;ic<NCoefficients;ic++){
		answer+=Mtot[ic]*A[ic];
	}
	return answer;
}
