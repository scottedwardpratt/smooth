#include <complex>
#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
#include "msu_smoothutils/log.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

void CMCMC::EvaluateTrace(){
	unsigned int itrace,ipar,jpar,ntrace=trace.size();
	Eigen::VectorXd thetabar;
	Eigen::MatrixXd sigma2;
	Eigen::VectorXcd evals;
	Eigen::MatrixXcd evecs;
	thetabar.resize(NPars);
	sigma2.resize(NPars,NPars);
	evecs.resize(NPars,NPars);
	thetabar.setZero();
	sigma2.setZero();

	
	
	
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			thetabar(ipar)+=trace[itrace].Theta[ipar];
			for(jpar=0;jpar<NPars;jpar++){
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
		for(jpar=0;jpar<NPars;jpar++){
			sigma2(ipar,jpar)=sigma2(ipar,jpar)-thetabar(ipar)*thetabar(jpar);
			//CLog::Info(to_string(sigma2(ipar,jpar))+" ");
		}
		//CLog::Info("\n");
	}
	cout << sigma2 << endl;
	
	Eigen::EigenSolver<Eigen::MatrixXd> esolver(sigma2);
	evals=esolver.eigenvalues();
	evecs=esolver.eigenvectors();
	Eigen::MatrixXcd sigma2test;
	Eigen::VectorXcd th,r;
	sigma2test.resize(NPars,NPars);
	th.resize(NPars);
	r.resize(NPars);
	sigma2test.setZero();
	unsigned int itest,ntest=10000000;
	vector<double> evalnorm;
	evalnorm.resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
		evalnorm[ipar]=sqrt(fabs(real(evals(ipar))));
	}
	
	printf("eigenvectors:\n");
	cout << evecs << endl;
	
	for(itest=0;itest<ntest;itest++){
		for(ipar=0;ipar<NPars;ipar++)
			r(ipar)=evalnorm[ipar]*randy->ran_gauss();
		th=evecs*r;
		for(ipar=0;ipar<NPars;ipar++){
			for(jpar=0;jpar<NPars;jpar++){
				sigma2test(ipar,jpar)+=th(ipar)*th(jpar);
			}
		}
	}
	sigma2test=sigma2test/double(ntest);
	
	cout << sigma2test << endl;
	
}

void CMCMC::OptimizeSteps(){
	unsigned int itrace,ipar,jpar,ntrace=trace.size();
	Eigen::VectorXd thetabar;
	thetabar.resize(NPars);
	thetabar.setZero();
	dThetadTheta.setZero();

	
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			thetabar(ipar)+=trace[itrace].Theta[ipar];
			for(jpar=0;jpar<NPars;jpar++){
				dThetadTheta(ipar,jpar)+=trace[itrace].Theta[ipar]*trace[itrace].Theta[jpar];
			}
		}
	}
	dThetadTheta=dThetadTheta/double(ntrace);
	thetabar=thetabar/double(ntrace);
	
	CLog::Info("<theta>=");
	for(ipar=0;ipar<NPars;ipar++){
		CLog::Info(to_string(thetabar(ipar))+" ");
	}
	CLog::Info("\n");
	
	CLog::Info("Sigma^2=\n");
	for(ipar=0;ipar<NPars;ipar++){
		for(jpar=0;jpar<NPars;jpar++){
			dThetadTheta(ipar,jpar)=dThetadTheta(ipar,jpar)-thetabar(ipar)*thetabar(jpar);
			//CLog::Info(to_string(sigma2(ipar,jpar))+" ");
		}
		//CLog::Info("\n");
	}
	cout << dThetadTheta << endl;
	
	Eigen::EigenSolver<Eigen::MatrixXd> esolver(dThetadTheta);
	dTdTEigenVals=esolver.eigenvalues();
	dTdTEigenVecs=esolver.eigenvectors();
	
	
	Eigen::MatrixXcd sigma2test;
	Eigen::VectorXcd th,r;
	sigma2test.resize(NPars,NPars);
	th.resize(NPars);
	r.resize(NPars);

	sigma2test.setZero();
	unsigned int itest,ntest=10000000;
	vector<double> evalnorm;
	evalnorm.resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
		evalnorm[ipar]=sqrt(fabs(real(dTdTEigenVals(ipar))));
	}
	
	printf("eigenvectors:\n");
	cout << dTdTEigenVecs << endl;
	
	for(itest=0;itest<ntest;itest++){
		for(ipar=0;ipar<NPars;ipar++)
			r(ipar)=evalnorm[ipar]*randy->ran_gauss();
		th=dTdTEigenVecs*r;
		for(ipar=0;ipar<NPars;ipar++){
			for(jpar=0;jpar<NPars;jpar++){
				sigma2test(ipar,jpar)+=th(ipar)*th(jpar);
			}
		}
	}
	sigma2test=sigma2test/double(ntest);
	
	cout << sigma2test << endl;
	
}