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

void CSmooth::CalcYDYDTheta(vector<double> &A,double LAMBDA,vector<double> &theta,double &Y,vector<double> &dYdTheta){
	unsigned int ir,ic,ipar,n,ndiff,oldipar,npar;
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
		iparvec.clear();
		nparvec.clear();
		ndiff=0;
		
		
		printf("------ ic=%u --------\n",ic);
		for(ir=0;ir<rank[ic];ir++){
			printf("%2u ",IPar[ic][ir]);
		}
		printf("\n");

		oldipar = 99999;

		for(ir = 0; ir < rank[ic]; ir++) {
			ipar = IPar[ic][ir];
			if(ipar != oldipar){
				iparvec.push_back(ipar);
				nparvec.push_back(1);
				//iparvec[ndiff] += 1;
				ndiff += 1;
				oldipar = ipar;
			}
			else {
				nparvec[ndiff-1]+=1 ;
				printf("check, ndiff=%u, nparvec=%u\n",ndiff,nparvec[ndiff-1]);
				oldipar=ipar;
			}
		}
		for(n=0;n<ndiff;n++){
			printf("%2u,%2u   ",iparvec[n],nparvec[n]);
		}
		printf("\n");
			
		printf("ndiff=%u, prefactor=%g, iparvec size=%lu, nparvec size=%lu\n",
		ndiff,prefactor,iparvec.size(),nparvec.size());
		
		term=prefactor;
		for(n = 0;n<ndiff;n++){
			ipar = iparvec[n];
			npar=nparvec[n];
			term*=pow(theta[ipar]/LAMBDA,npar);
		}
		
		for(n=0;n<ndiff;n++){
			ipar = iparvec[n];
			npar=nparvec[n];
			dYdTheta[ipar]+=term*npar/theta[ipar];
		}
		
	}	
	for(unsigned int i = 0; i < dYdTheta.size(); i++) {
		dYdTheta[i] *= rfactor;
	}
}

