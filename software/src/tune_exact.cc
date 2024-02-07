#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;


void CSmoothEmulator::TuneExact(){
	unsigned int NCoefficients=smooth->NCoefficients;
	unsigned int itrain,ic,a,b,c;
	CalcMForTraining();
	Eigen::VectorXd alpha,gamma,YTrain,delta;
	Eigen::MatrixXd beta,C;
	alpha.resize(NCoefficients);
	gamma.resize(NCoefficients);
	beta.resize(NTrainingPts,NCoefficients);
	C.resize(NTrainingPts,NTrainingPts);
	delta.resize(NTrainingPts);
	alpha.setZero();
	beta.setZero();
	gamma.setZero();
	C.setZero();
	delta.setZero();
	Abest.resize(NCoefficients);
	double betaadotbetab,betaadotbetac,betabdotbetac;
	YTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTrain[itrain]=smoothmaster->traininginfo->YTrain[iY][itrain];
	}
	alpha=Minv*YTrain;
	for(a=0;a<NTrainingPts;a++){
		for(b=NTrainingPts;b<NCoefficients;b++){
			for(itrain=0;itrain<NCoefficients;itrain++){
				beta(a,b)+=Minv(a,itrain)*Mtot[itrain][b];
			}
		}
	}
	// Will solve for gamma where delta=C*gamma, but first find C
	for(b=0;b<NTrainingPts;b++){
		for(c=0;c<NTrainingPts;c++){
			betabdotbetac=0.0;
			for(ic=NTrainingPts;ic<NCoefficients;ic++){
				betabdotbetac+=beta(b,ic)*beta(c,ic);
			}
			C(b,c)+=betabdotbetac;
			for(a=0;a<NTrainingPts;a++){
				betaadotbetab=betaadotbetac=0.0;
				for(ic=NTrainingPts;ic<NCoefficients;ic++){
					betaadotbetab+=beta(a,ic)*beta(b,ic);
					betaadotbetac+=beta(a,ic)*beta(c,ic);
				}
				C(b,c)+=betaadotbetab*betaadotbetac;
			}
		}
	}
	for(b=0;b<NTrainingPts;b++){
		delta(b)=0.0;
		for(a=0;a<NTrainingPts;a++){
			betaadotbetab=0.0;
			for(ic=NTrainingPts;ic<NCoefficients;ic++){
				betaadotbetab+=beta(a,ic)*beta(b,ic);
			}
			delta(b)+=betaadotbetab*alpha(a);
		}
	}
	gamma=C.colPivHouseholderQr().solve(delta);
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		Abest[ic]=0.0;
		for(a=0;a<NTrainingPts;a++){
			Abest[ic]+=gamma(a)*beta(a,ic);
		}
	}
	for(a=0;a<NTrainingPts;a++){
		Abest[a]=alpha(a);
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			Abest[a]-=beta(a,ic)*Abest[ic];
		}
	}
}
