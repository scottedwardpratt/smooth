#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;

CSmoothEmulator::TuneExact(){
	CalcMForTraining();
	observable_name=observable_name_set;
	NTrainingPts=smoothmaster->traininginfo->NTrainingPts;

}




void CSmoothEmulator::SetThetaTrain(){
	ThetaTrain.resize(NTrainingPts);
	for(unsigned int itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(unsigned int ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=smoothmaster->traininginfo->modelpars[itrain]->Theta[ipar];
		}
	}
}

void CSmoothEmulator::Init(){
	SigmaA=SigmaA0;
	NSigmaA=0;
	SigmaAbar=0.0;
	FirstTune=true;
	MCStepSize=MCStepSize/double(NPars*NPars);
	MCSigmaAStepSize=MCSigmaAStepSize/double(NPars*NPars);

	ASample.resize(NASample);
	for(unsigned int isample=0;isample<NASample;isample++){
		ASample[isample].resize(smooth->NCoefficients);
		//SetA_RanGauss(SigmaA,ASample[isample]);
		SetA_Zero(ASample[isample]);
	}

	A.resize(smooth->NCoefficients);
	SetA_Zero(A);
	ATrial.resize(smooth->NCoefficients);
	SetA_Zero(ATrial);
	M.resize(NTrainingPts,NTrainingPts);
	Minv.resize(NTrainingPts,NTrainingPts);
	Mtot.resize(NTrainingPts);
	for(unsigned int it=0;it<NTrainingPts;it++){
		Mtot[it].resize(smooth->NCoefficients);
		for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
			Mtot[it][ic]=0.0;
		}
	}
}

void CSmoothEmulator::Tune(){
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
	double betadotbetab,betaadotbetac,betabdotbetac;
	YTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPtsitrain++){
		YTrain[itrain]=smoothmaster->traininginfo->YTrain[iY][itrain];
	}
	alpha=Minv*YTrain;
	for(a=0;a<NTrain;a++){
		for(b=NTrain;b<NCoefficients;b++){
			for(itrain=0;itrain<NCoefficients;itrain++){
				beta(a,b)+=Minv(a,itrain)*Mtot(itrain,b);
			}
		}
	}
	// Will solve for gamma where delta=C*gamma, but first find C
	for(b=0;b<NTrainingPts;b++){
		for(c=0;c<NTrainingPts;c++){
			betabdotbetac=0.0;
			for(ic=NTrainingPts;ic<NCoefficienta;ic++){
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
			delta(b)+=betadotbetab*alpha(a);
		}
	}
	gamma=C.colPivHouseholderQr().solve(delta);
	for(ic=NTrainingPts;ic<Coefficients;ic++){
		Abest(ic)=0.0;
		for(a=0;a<NTrainingPts;a++){
			Abest(ic)+=gamma(a)*beta(a,ic);
		}
	}
	for(a=0;a<NTrainingPts;a++){
		Abest(a)=alpha(a);
		for(ic=NTrainingPts;ic<NCoefficents;ic++){
			Abest(a)-=beta(a,ic)*Abest(ic);
		}
	}
}
