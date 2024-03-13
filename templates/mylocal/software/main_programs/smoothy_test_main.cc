#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

int main(){
	unsigned int maxrank=4;
	unsigned int npars=14;
	unsigned int ic,itheta;
	//unsigned int idtheta;
	CSmooth *smoothy=new CSmooth(npars,maxrank);
	unsigned int Ncoefficients=smoothy->NCoefficients;
	double lambda=2.5,Y;
	//double Y0,Y1,Y2,dtheta=0.01,error=0.0;
	vector<double> A(Ncoefficients),dYdTheta(Ncoefficients);
	vector<double> theta(npars),theta1(npars),theta2(npars);
	Crandy *randy=new Crandy(1234);
	
	unsigned int itrial,ntrial=100000;
	
	for(itrial=0;itrial<ntrial;itrial++){
		for(ic=0;ic<Ncoefficients;ic++){
			A[ic]=100.0*randy->ran();
		}
	
		for(itheta=0;itheta<npars;itheta++){
			theta[itheta]=-1.0+2.0*randy->ran();
		}
		//Y=smoothy->CalcY(A,lambda,theta);
		smoothy->CalcYDYDTheta(A,lambda,theta,Y,dYdTheta);
		/*
		Y0=Y;
		smoothy->CalcYDYDTheta(A,lambda,theta,Y,dYdTheta);
		
		error=fabs(Y0-Y);
		for(idtheta=0;idtheta<npars;idtheta++){
			for(itheta=0;itheta<npars;itheta++){
				theta1[itheta]=theta[itheta];
				theta2[itheta]=theta[itheta];
				if(itheta==idtheta){
					theta1[itheta]-=0.5*dtheta;
					theta2[itheta]+=0.5*dtheta;
				}
			}
			Y2=smoothy->CalcY(A,lambda,theta2);
			Y1=smoothy->CalcY(A,lambda,theta1);
			//printf("%g =? %g\n",dYdTheta[idtheta],(Y2-Y1)/dtheta);
			error+=fabs(dYdTheta[idtheta]-(Y2-Y1)/dtheta);
		}
		*/
		if((itrial+1)%(ntrial/10)==0){
			printf("finished %g percent\n",100.0*(itrial+1)/double(ntrial));
			//error=0.0;
		}
		
		
	}
	return 0;
}
