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

void CSmooth::CalcDYDTheta(vector<double> &A,double LAMBDA,vector<double> &theta,double &Y,vector<double> &dYdTheta){
	unsigned int ir,ic;
	dYdTheta.resize(theta.size(),0.0);
	double rfactor = GetRFactor(LAMBDA, theta);
	double term,prefactor;

	vector<unsigned int> nparvec;
	vector<unsigned int> iparvec;
	nparvec.resize(0);
	iparvec.resize(0);
	
	Y=0.0;
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
		}
		Y+=term;
	}
	Y*=rfactor;
	

	for(unsigned int ic = 0; ic < NCoefficients; ic++) {
		prefactor = A[ic] * sqrt(double(dupfactor[ic]) / double(factorial[rank[ic]]));

		unsigned int oldipar = -1;
		unsigned int ndiff = 0;

		for(ir = 0; ir < rank[ic]; ir++) {
			unsigned ipar = IPar[ic][ir];
			if(ipar != oldipar){
				iparvec.push_back(ipar);
				iparvec[ndiff] += 1;
				ndiff += 1;
				oldipar = ipar;
			}
			else {
				nparvec[ndiff]+=1 ;
			}

			for(unsigned int n = 0;n<ndiff;n++){
				ipar = iparvec[n];
				term = prefactor;
				for(unsigned int nprime = 0;nprime< ndiff;nprime++){

					unsigned int iparprime = iparvec[nprime];
					unsigned int nprime_val = nparvec[iparprime];
					if(iparprime != ipar){

						term *= pow(theta[iparprime],nprime_val);
					}
					else{
						term *= nprime_val * pow(theta[ipar]/LAMBDA,nprime_val-1);
					}
				}
				dYdTheta[ipar] += term/LAMBDA;
			}
		}
	}
	for(unsigned int i = 0; i < dYdTheta.size(); i++) {
		dYdTheta[i] *= rfactor;
	}
}

