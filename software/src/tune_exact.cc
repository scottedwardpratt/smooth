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
	H6.resize(NTrainingPts);
	H8.resize(NTrainingPts);
	Psi.resize(NTrainingPts,NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		BetaDotBeta[a].resize(NTrainingPts);
		Mdotbeta[a].resize(NTrainingPts);
		H6[a].resize(NTrainingPts);
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
	AExact.resize(NCoefficients);
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
		AExact[ic]=0.0;
		for(a=0;a<NTrainingPts;a++){
			AExact[ic]+=gamma(a)*beta(a,ic);
		}
	}

	for(a=0;a<NTrainingPts;a++){
		AExact[a]=alpha(a);
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			AExact[a]-=beta(a,ic)*AExact[ic];
		}
	}
	
	// Now calculate arrays used for calculating uncertainty
	
	Eigen::MatrixXd Psi,D;
	Psi.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	D.setZero();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			if(a==b)
				D(a,b)=1.0;
			D(a,b)+=BetaDotBeta[a][b];
		}
	}

	Psi=-D.inverse();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			H6[a][b]=0.0;
			H8[a][b]=0.0;
			for(aprime=0;aprime<NTrainingPts;aprime++){
				H6[a][b]+=BetaDotBeta[a][aprime]*Psi(aprime,b);
				for(bprime=0;bprime<NTrainingPts;bprime++){
					H8[a][b]+=BetaDotBeta[a][aprime]*Psi(aprime,bprime)*BetaDotBeta[bprime][b];
				}
			}
		}
	}	

}

void CSmoothEmulator::GetExactUncertainty(vector<double> &Theta_s,double &sigma){
	double sigma2;
	unsigned int ic,a,b,NCoefficients=smooth->NCoefficients;
	Eigen::VectorXd Mtot_s,M_s,Mtot_sDotBeta;
	M_s.resize(NTrainingPts);
	Mtot_s.resize(NCoefficients);
	M_s.setZero();
	Mtot_s.setZero(NCoefficients);
	Mtot_sDotBeta.resize(NTrainingPts);
	Mtot_sDotBeta.setZero();
	
	for(a=0;a<NTrainingPts;a++){
		M_s[a]=smooth->GetM(a,LAMBDA,Theta_s);
	}
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		Mtot_s[ic]=smooth->GetM(ic,LAMBDA,Theta_s);
	}
	
	for(a=0;a<NTrainingPts;a++){
		Mtot_sDotBeta(a)=0.0;
		for(ic=NTrainingPts;ic<NCoefficients;ic++){
			Mtot_sDotBeta(a)+=Mtot_s(ic)*beta(a,ic);
		}
	}

	sigma2=0.0;
	// First term
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		sigma2+=Mtot_s(ic)*Mtot_s(ic);
	}
	// Second  & third terms
	for(a=0;a<NTrainingPts;a++){
		sigma2-=2*M_s(a)*Mtot_sDotBeta(a);
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
			sigma2-=2*M_s(a)*H6[a][b]*Mtot_sDotBeta(b);
		}
	}
	/// 8th term
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			sigma2+=M_s(a)*H8[a][b]*M_s(b);
		}
	}
	printf("---- sigma2=%g\n",sigma2);
	if(sigma2<1.0E-8){
		CLog::Info("Inside CSmoothEmulator::GetExactUncertainty, sigma^2 is less than zero = "+to_string(sigma2)+"\n");
	}
	sigma=sqrt(fabs(sigma2));	

}


