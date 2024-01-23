#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

void CMCMC::EvaluateTrace(){
	unsigned int itrace,ipar,jpar,ntrace=trace.size();
	vector<double> thetabar(ntrace,0.0);
	vector<vector<double>> sigma2;
	sigma2.resize(ntrace);
	for(ipar=0;ipar<NPars;ipar++){
		sigma2[ipar].resize(NPars);
		for(jpar=0;jpar<NPars;jpar++)
			sigma2[ipar][jpar]=0.0;
	}
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			thetabar[ipar]+=trace[itrace].Theta[ipar];
			for(jpar=0;jpar<=ipar;jpar++){
				sigma2[ipar][jpar]+=trace[itrace].Theta[ipar]*trace[itrace].Theta[jpar];
			}
		}
	}
	for(ipar=0;ipar<NPars;ipar++){
		thetabar[ipar]=thetabar[ipar]/double(ntrace);
		for(jpar=0;jpar<=ipar;jpar++){
			sigma2[ipar][jpar]=sigma2[ipar][jpar]/double(ntrace);
		}
	}
	for(ipar=0;ipar<NPars;ipar++){
		printf("<theta[%u]>=%g\n",ipar,thetabar[ipar]);
	}
	printf("Sigma^2=\n");
	for(ipar=0;ipar<NPars;ipar++){
		for(jpar=0;jpar<=ipar;jpar++){
			sigma2[ipar][jpar]=sigma2[ipar][jpar]-thetabar[ipar]*thetabar[jpar];
			printf("%9.5f ",sigma2[ipar][jpar]);
		}
		printf("\n");
	}
	
}