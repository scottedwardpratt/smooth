#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::SetSigmaA(double sigmaAset){
	SigmaA=sigmaAset;
	FixSigmaA=true;
	Bcalculated=false;
}

void CSmoothEmulator::CalcSigmaA(){
	int a,b;
	/*
	double sigmaA2=0.0;
	for(a=0;a<int(NTrainingPts);a++){
	for(b=0;b<int(NTrainingPts);b++){
	sigmaA2+=smoothmaster->traininginfo->YTrain[iY][a]
	*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
	}
	}
	sigmaA2=sigmaA2/double(NTrainingPts);
	SigmaA=sqrt(sigmaA2);
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	*/
	
	/*
	double TrB,SumB;
	TrB=SumB=0.0;
	for(a=0;a<int(NTrainingPts);a++){
	for(b=0;b<int(NTrainingPts);b++){
	SumB+=B(a,b);
	}
	TrB+=B(a,a);
	}
	double ysquared,Ya,ybar;
	ysquared=ybar=0.0;
	for(a=0;a<int(NTrainingPts);a++){
	Ya=smoothmaster->traininginfo->YTrain[iY][a];
	ysquared+=Ya*Ya;
	ybar+=Ya;
	}
	ybar=ybar/double(NTrainingPts);
	SigmaA=(ysquared-ybar*ybar)/(TrB-SumB/double(NTrainingPts*NTrainingPts));
	SigmaA=sqrt(SigmaA);
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	*/
	
	
	double TrB,SumB;
	TrB=SumB=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		for(b=0;b<int(NTrainingPts);b++){
			SumB+=B(a,b);
		}
		TrB+=B(a,a);
	}
	double ysquared,Ya,ybar;
	ysquared=ybar=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		Ya=smoothmaster->traininginfo->YTrain[iY][a];
		ysquared+=Ya*Ya;
		ybar+=Ya;
	}
	ybar=ybar/double(NTrainingPts);
	SigmaA=ysquared/TrB;
	SigmaA=sqrt(SigmaA);
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	

}

void CSmoothEmulator::GetSigmaA123(double &sig1,double &sig2,double &sig3){
	int a,b;
	
	double sigmaA2=0.0;
	for(a=0;a<int(NTrainingPts);a++){
	for(b=0;b<int(NTrainingPts);b++){
	sigmaA2+=smoothmaster->traininginfo->YTrain[iY][a]
	*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
	}
	}
	sigmaA2=sigmaA2/double(NTrainingPts);
	SigmaA=sqrt(sigmaA2);
	sig1=SigmaA;
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	
	double TrB,SumB;
	TrB=SumB=0.0;
	for(a=0;a<int(NTrainingPts);a++){
	for(b=0;b<int(NTrainingPts);b++){
	SumB+=B(a,b);
	}
	TrB+=B(a,a);
	}
	double ysquared,Ya,ybar;
	ysquared=ybar=0.0;
	for(a=0;a<int(NTrainingPts);a++){
	Ya=smoothmaster->traininginfo->YTrain[iY][a];
	ysquared+=Ya*Ya;
	ybar+=Ya;
	}
	ybar=ybar/double(NTrainingPts);
	SigmaA=(ysquared-ybar*ybar)/(TrB-SumB/double(NTrainingPts*NTrainingPts));
	SigmaA=sqrt(SigmaA);
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	sig2=SigmaA;
	
	TrB=SumB=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		for(b=0;b<int(NTrainingPts);b++){
			SumB+=B(a,b);
		}
		TrB+=B(a,a);
	}

	ysquared=ybar=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		Ya=smoothmaster->traininginfo->YTrain[iY][a];
		ysquared+=Ya*Ya;
		ybar+=Ya;
	}
	ybar=ybar/double(NTrainingPts);
	SigmaA=ysquared/TrB;
	SigmaA=sqrt(SigmaA);
	sig3=SigmaA;
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	

}

void CSmoothEmulator::CalcExactLogP(){
	double exparg=0.0;
	for(int a=0;a<int(NTrainingPts);a++){
		for(int b=0;b<int(NTrainingPts);b++){
			exparg-=0.5*smoothmaster->traininginfo->YTrain[iY][a]*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	//printf("exparg^1/2=%g\n",sqrt(fabs(exparg)));
	exparg=exparg/(SigmaA*SigmaA);
	double detB=B.determinant();
	logP=0.5*log(fabs(detB))-NTrainingPts*log(SigmaA)+exparg;
	
}

void CSmoothEmulator::CalcSigmaLambda(){
	int a,alpha,beta;
	double B,ybar=0.0;
	Eigen::MatrixXd TT,TTinv;
	Eigen::VectorXd M,TYbar,x;
	TT.resize(NPars,NPars);
	TTinv.resize(NPars,NPars);
	M.resize(NPars);
	TYbar.resize(NPars);
	TT.setZero();
	for(a=0;a<int(NTrainingPts);a++){
		ybar+=smoothmaster->traininginfo->YTrain[iY][a];
		for(alpha=0;alpha<int(NPars);alpha++){
			for(beta=0;beta<int(NPars);beta++){
				TT(alpha,beta)+=smoothmaster->traininginfo->modelpars[a]->Theta[alpha]*smoothmaster->traininginfo->modelpars[a]->Theta[beta];
			}
		}
	}
	TT=TT/((1.0/double(NTrainingPts))-1.0);
	ybar=ybar/double(NTrainingPts);
	TTinv=TT.inverse();
	for(alpha=0;alpha<int(NPars);alpha++){
		TYbar[alpha]=0.0;
		for(a=0;a<int(NTrainingPts);a++)
			TYbar[alpha]+=ybar*smoothmaster->traininginfo->modelpars[a]->Theta[alpha];
	}
	M=TTinv*TYbar;
	B=ybar;
	for(a=0;a<int(NTrainingPts);a++){
		for(alpha=0;alpha<int(NPars);alpha++){
			B-=(1.0/double(NTrainingPts))*smoothmaster->traininginfo->modelpars[a]->Theta[alpha]*M(alpha);
		}
	}
	printf("b=%g, M=(",B);
	for(alpha=0;alpha<int(NTrainingPts);alpha++){
		printf("%g ",M(alpha));
	}
	printf(")\n");
	
	
	
	
	
}
