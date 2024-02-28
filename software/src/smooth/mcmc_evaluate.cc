#include <complex>
#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
#include "msu_smoothutils/log.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

void CMCMC::EvaluateTrace(){
	unsigned int itrace,ipar,jpar,ntrace=trace.size();
	vector<double> thetabar;
	string SigmaString;
	char cc[CLog::CHARLENGTH];
	Eigen::MatrixXd sigma2;
	Eigen::VectorXcd evals;
	Eigen::MatrixXcd evecs;
	thetabar.resize(NPars);
	sigma2.resize(NPars,NPars);
	evecs.resize(NPars,NPars);
	sigma2.setZero();
	for(ipar=0;ipar<NPars;ipar++){
		thetabar[ipar]=0.0;
	}
	
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			thetabar[ipar]+=trace[itrace].Theta[ipar];
			for(jpar=0;jpar<NPars;jpar++){
				sigma2(ipar,jpar)+=trace[itrace].Theta[ipar]*trace[itrace].Theta[jpar];
			}
		}
	}
	sigma2=sigma2/double(ntrace);
	for(ipar=0;ipar<NPars;ipar++){
		thetabar[ipar]=thetabar[ipar]/double(ntrace);
	}
	
	CModelParameters modpars;
	modpars.priorinfo=master->priorinfo;
	modpars.SetTheta(thetabar);
	modpars.Print();
	
	CLog::Info("Sigma^2=\n");
	for(ipar=0;ipar<NPars;ipar++){
		SigmaString.clear();
		for(jpar=0;jpar<NPars;jpar++){
			sigma2(ipar,jpar)=sigma2(ipar,jpar)-thetabar[ipar]*thetabar[jpar];
			snprintf(cc,CLog::CHARLENGTH,"%12.5e ",sigma2(ipar,jpar));
			SigmaString=SigmaString+cc;
		}
		SigmaString=SigmaString+"\n";
		CLog::Info(SigmaString);
	}
	
	Eigen::EigenSolver<Eigen::MatrixXd> esolver(sigma2);
	evals=esolver.eigenvalues();
	evecs=esolver.eigenvectors();
	Eigen::VectorXcd th,r;
	th.resize(NPars);
	r.resize(NPars);
	vector<double> evalnorm;
	evalnorm.resize(NPars);
	CLog::Info("EigenValues for Sigma:\n");
	for(ipar=0;ipar<NPars;ipar++){
		evalnorm[ipar]=sqrt(fabs(real(evals(ipar))));
		CLog::Info(to_string(evalnorm[ipar])+" \n");
	}
	
	printf("eigenvectors:\n");
	for(ipar=0;ipar<NPars;ipar++){
		SigmaString.clear();
		for(jpar=0;jpar<NPars;jpar++){
			snprintf(cc,CLog::CHARLENGTH,"%12.5e ",real(evecs(ipar,jpar)));
			SigmaString+=cc;
		}
		SigmaString+="\n";
		CLog::Info(SigmaString);
	}	
}