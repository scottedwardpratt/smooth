#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
#include "msu_commonutils/log.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

void CMCMC::EvaluateTrace(){
	unsigned int itrace,ipar,jpar,ntrace=trace.size();
	Eigen::VectorXd thetabar;
	Eigen::MatrixXd sigma2;
	thetabar.resize(ntrace);
	sigma2.resize(ntrace,ntrace);
	thetabar.setZero();
	sigma2.setZero();
	
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			thetabar(ipar)+=trace[itrace].Theta[ipar];
			for(jpar=0;jpar<=ipar;jpar++){
				sigma2(ipar,jpar)+=trace[itrace].Theta[ipar]*trace[itrace].Theta[jpar];
			}
		}
	}
	sigma2=sigma2/double(ntrace);
	thetabar=thetabar/double(ntrace);
	
	CLog::Info("<theta>=");
	for(ipar=0;ipar<NPars;ipar++){
		CLog::Info(to_string(thetabar(ipar))+" ");
	}
	CLog::Info("\n");
	
	CLog::Info("Sigma^2=\n");
	for(ipar=0;ipar<NPars;ipar++){
		for(jpar=0;jpar<=NPars;jpar++){
			sigma2(ipar,jpar)=sigma2(ipar,jpar)-thetabar(ipar)*thetabar(jpar);
			CLog::Info(to_string(sigma2(ipar,jpar))+" ");
		}
		CLog::Info("\n");
	}
	
	
	
}