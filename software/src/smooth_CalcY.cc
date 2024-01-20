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
	unsigned int ir,ic,ipar,iparprime,n,nprime,nparprime,ndiff,oldipar,npar;
	dYdTheta.resize(theta.size(),0.0);
	double rfactor=GetRFactor(LAMBDA, theta);
	double term,prefactor;

	vector<unsigned int> nparvec;
	vector<unsigned int> iparvec;
	nparvec.resize(0);
	iparvec.resize(0);
	dYdTheta.resize(NPars);
	for(ipar=0;ipar<NPars;ipar++)
		dYdTheta[ipar]=0.0;
	Y=0.0;
	
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		if(A[ic]!=A[ic]){
			CLog::Fatal("Disaster in CSmooth::CalcYDYDTheta, A!=A\n");
		}
		if(term!=term){
			CLog::Fatal("Disaster in CSmooth::CalcYDYDTheta, term!=term\n");
		}
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
		}
		Y+=term;
	}
	Y*=rfactor;

	double termtest,Ytest=0.0;
	for(unsigned int ic=0;ic<NCoefficients;ic++){
		prefactor= A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		iparvec.clear();
		nparvec.clear();
		ndiff=0;
		oldipar=99999;
		for(ir=0;ir<rank[ic];ir++){
			ipar=IPar[ic][ir];
			if(ipar!=oldipar){
				iparvec.push_back(ipar);
				nparvec.push_back(1);
				ndiff+=1;
				oldipar=ipar;
			}
			else{
				nparvec[ndiff-1]+=1 ;
			}
		}
		
		termtest=prefactor;
		for(n=0;n<ndiff;n++){
			ipar=iparvec[n];
			npar=nparvec[n];
			termtest*=pow(theta[ipar]/LAMBDA,npar);
			term=prefactor;
			for(nprime=0;nprime<ndiff;nprime++){
				iparprime=iparvec[nprime];
				nparprime=nparvec[nprime];
				if(nprime!=n)
					term*=pow(theta[iparprime]/LAMBDA,nparprime);
				else
					term*=(double(npar)/LAMBDA)*pow(theta[ipar]/LAMBDA,npar-1);
			}
			dYdTheta[ipar]+=term;
			if(dYdTheta[ipar]!=dYdTheta[ipar]){
				printf("ipar=%u, npar=%u, term=%g, theta=%g\n",ipar,npar,term,theta[ipar]);
				CLog::Fatal("Disaster in CSmooth::CalcYDYDTheta\n");
			}
		}
		Ytest+=termtest;
	}	
	printf("Y=%g, Ytest=%g\n",Y,Ytest*rfactor);
	for(ipar=0;ipar<dYdTheta.size();ipar++) {
		dYdTheta[ipar]*=rfactor;
	}
}

