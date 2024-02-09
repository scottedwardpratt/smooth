#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;


void CSmoothEmulator::TuneExact(){
	unsigned int NCoefficients=smooth->NCoefficients;
	unsigned int itrain,ic,a,b,c,aprime,bprime;
	CalcMForTraining();
	BetaDotBeta.resize(NTrainingPts);
	Mdotbeta.resize(NTrainingPts);
	Psi.resize(NTrainingPts,NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		BetaDotBeta[a].resize(NTrainingPts);
		Mdotbeta[a].resize(NTrainingPts);
		H8[a].resize(NTrainingPts);
	}
	Eigen::VectorXd alpha,gamma,YTrain,delta;
	Eigen::MatrixXd C;
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
	YTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTrain[itrain]=smoothmaster->traininginfo->YTrain[iY][itrain];
	}
	alpha=Minv*YTrain;
	for(a=0;a<NTrainingPts;a++){
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			beta(a,ic)=0.0;
			for(b=0;b<NTrainingPts;b++){
				beta(a,ic)+=Minv(a,b)*Mtot[b][ic];
			}
		}
	}
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			BetaDotBeta[a][b]=0.0;
			Mdotbeta[a][b]=0.0;
			for(ic=NTrainingPts;ic<NCoefficients;ic++){
				BetaDotBeta[a][b]+=beta(a,ic)*beta(b,ic);
				Mdotbeta[a][b]+=Mtot[a][ic]*beta(b,ic);
			}
		}
	}
	// Will solve for gamma where delta=C*gamma, but first find C
	for(b=0;b<NTrainingPts;b++){
		for(c=0;c<NTrainingPts;c++){
			C(b,c)+=BetaDotBeta[b][c];
			for(a=0;a<NTrainingPts;a++){
				C(b,c)+=BetaDotBeta[a][b]*BetaDotBeta[a][c];
			}
		}
	}
	
	for(b=0;b<NTrainingPts;b++){
		delta(b)=0.0;
		for(a=0;a<NTrainingPts;a++){
			delta(b)+=BetaDotBeta[b][b]*alpha(a);
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
	
	// Now calculate arrays used for calculating uncertainty
	
	Eigen::MatrixXd Psi,D;
	Psi.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	D.setZero();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;a++){
			if(a==b)
				D(a,b)=1.0;
			D(a,b)+=BetaDotBeta[a][b];
		}
	}
	Psi=-D.inverse();
	for(a=0;a<NTrainingPts;a++){
		for(b=a;b<NTrainingPts;b++){
			H8[a][b]=0.0;
			for(aprime=0;aprime<NTrainingPts;aprime++){
				for(bprime=0;bprime<NTrainingPts;bprime++){
					H8[a][b]+=BetaDotBeta[a][aprime]*Psi(aprime,bprime)*BetaDotBeta[bprime][b];
				}
			}
		}
		if(a!=b)
			H8[b][a]=H8[a][b];
	}	
	
}

void CSmoothEmulator::GetExactUncertainty(vector<double> &Theta_s,double &sigma){
	double sigma2;
	unsigned int ic,a,b,NCoefficients=smooth->NCoefficients;
	Eigen::VectorXd Mtot_s,M_s,M_sDotBeta,Mtot_sDotBeta,M_sDotBetaDotBeta;
	M_s.resize(NTrainingPts);
	M_sDotBeta.resize(NCoefficients);
	Mtot_s.resize(NCoefficients);
	M_s.setZero();
	Mtot_s.setZero(NCoefficients);
	Mtot_sDotBeta.resize(NTrainingPts);
	Mtot_sDotBeta.setZero();
	M_sDotBetaDotBeta.resize(NTrainingPts);
	
	for(ic=0;ic<NCoefficients;ic++){
		Mtot_s[ic]=smooth->GetM(ic,LAMBDA,Theta_s);
		if(ic<NTrainingPts)
			M_s[ic]=Mtot_s[ic];
	}
	
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		M_sDotBeta(ic)=0.0;
		for(a=0;a<NTrainingPts;a++){
			M_sDotBeta(ic)+=M_s[a]*beta(a,ic);
		}
	}
	
	for(a=0;a<NTrainingPts;a++){
		Mtot_sDotBeta(a)=0.0;
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			Mtot_sDotBeta(a)+=Mtot_s(ic)*beta(a,ic);
		}
	}
	
	M_sDotBetaDotBeta.setZero();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			M_sDotBetaDotBeta(a)+=M_s(b)*BetaDotBeta[a][b];
		}
	}

	sigma2=0.0;
	// First term
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		sigma2+=Mtot_s(ic)*Mtot_s(ic);
	}
	// Second  & third terms
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		sigma2-=2*M_sDotBeta(ic)*Mtot_s(ic);
	}
	// Fourth term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			sigma2+=M_s(a)*BetaDotBeta[a][b]*M_s(b);
		}
	}
	// Fifth term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			sigma2+=Mtot_sDotBeta(a)*Psi(a,b)*Mtot_sDotBeta(b);
		}
	}
	// Sixth & seventh terms
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			sigma2-=2*M_sDotBetaDotBeta(a)*Psi(a,b)*Mtot_sDotBeta(b);
		}
	}
	/// 8th term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			sigma2+=M_s(a)*H8[a][b]*M_s(b);
		}
	}
	
	if(sigma2<0.0){
		CLog::Fatal("Inside CSmoothEmulator::GetExactUncertainty, sigma^2 is less than zero = "+to_string(sigma2)+"\n");
	}
	sigma=sqrt(sigma2);	

}


