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
	
	double sigmaA2=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		for(b=0;b<int(NTrainingPts);b++){
			sigmaA2+=smoothmaster->traininginfo->YTrain[iY][a]
				*Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	sigmaA2=sigmaA2/double(NTrainingPts);
	SigmaA=sqrt(fabs(sigmaA2));
	
	//CLog::Info("SigmaA="+to_string(SigmaA)+"\n");
	
	
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
	*/
	

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
	exparg=exparg/(SigmaA*SigmaA);
	double detB=B.determinant();
	logP=-0.5*log(fabs(detB))-NTrainingPts*log(SigmaA)+exparg;
	if(isinf(logP)){
		printf("XXXXXXXXX detB=%g, Lambda=%g, SigmaA=%g, exparg=%g\n",detB,LAMBDA,SigmaA,exparg);
	}
	
}

void CSmoothEmulator::CalcSigmaLambdaAlt(double &LambdaGuess){
	int a,alpha,beta;
	double B,ybar;
	vector<vector<double>> dTheta;
	vector<double> dy;
	Eigen::MatrixXd TT,TTinv;
	Eigen::VectorXd M,Thetabar,TY;
	dy.resize(NTrainingPts);
	dTheta.resize(NTrainingPts);
	Thetabar.resize(NTrainingPts);
	for(a=0;a<int(NTrainingPts);a++){
		dTheta[a].resize(NPars);
		for(alpha=0;alpha<int(NPars);alpha++)
			dTheta[a][alpha]=0.0;
	}
	
	ybar=0.0;
	for(a=0;a<int(NTrainingPts);a++){
		dy[a]=smoothmaster->traininginfo->YTrain[iY][a];
		ybar+=dy[a];
	}
	ybar=ybar/double(NTrainingPts);
	for(a=0;a<int(NTrainingPts);a++)
		dy[a]-=ybar;
	
	for(alpha=0;alpha<int(NPars);alpha++){
		Thetabar[alpha]=0.0;
		for(a=0;a<int(NTrainingPts);a++){
			dTheta[a][alpha]=smoothmaster->traininginfo->modelpars[a]->Theta[alpha];
			Thetabar[alpha]+=dTheta[a][alpha];
		}
		Thetabar[alpha]=Thetabar[alpha]/double(NTrainingPts);
		for(a=0;a<int(NTrainingPts);a++){
			dTheta[a][alpha]-=Thetabar[alpha];
		}
	}
	
	TT.resize(NPars,NPars);
	TTinv.resize(NPars,NPars);
	M.resize(NPars);
	TY.resize(NPars);
	TT.setZero();
	
	
	for(alpha=0;alpha<int(NPars);alpha++){
		TY(alpha)=0.0;
		for(a=0;a<int(NTrainingPts);a++){
			TY(alpha)+=dy[a]*dTheta[a][alpha];
		}
		for(beta=0;beta<int(NPars);beta++){
			TT(alpha,beta)=0.0;
			for(a=0;a<int(NTrainingPts);a++){
				TT(alpha,beta)+=dTheta[a][alpha]*dTheta[a][beta];
			}
		}
	}

	TTinv=TT.inverse();
	M=TTinv*TY;
	
	B=ybar;
	for(alpha=0;alpha<int(NPars);alpha++){
		B-=Thetabar[alpha]*M(alpha);
	}
	
	double residual=0.0,yguess;
	for(a=0;a<int(NTrainingPts);a++){
		yguess=B;
		for(alpha=0;alpha<int(NPars);alpha++){
			yguess+=M(alpha)*smoothmaster->traininginfo->modelpars[a]->Theta[alpha];
		}
		residual+=pow(yguess-smoothmaster->traininginfo->YTrain[iY][a],2);
	}
	//printf("residual=%g\n",residual);
	LambdaGuess=(NTrainingPts-NPars)*pow(residual,-1.0/3.0);
	//printf("LambdaGuess=%g\n",LambdaGuess);	
	
}

void CSmoothEmulator::CalcSigmaLambda(){
	double LambdaMin=2.0;
	double bestLambda,dLambda=1.0,bestlogP,oldbestlogP,oldbestLambda;
	int nfail=0;
	LAMBDA=LambdaMin;  // minimum LAMBDA
	Bcalculated=false;
	CalcBTTrain();
	CalcSigmaA();
	CalcExactLogP();
	bestLambda=LAMBDA;
	bestlogP=logP;
	oldbestlogP=logP;
	oldbestLambda=LAMBDA;
	if(logP!=logP || isinf(logP)){
		nfail=100;
		CLog::Info("Note: LAMBDA set to "+to_string(LAMBDA)+" -- optimum value might be higher\n but determinant(B) cannot be calculated due to numerical accuracy problems.\n");
		logP=-200.0;
	}
	LAMBDA=bestLambda+dLambda;
	while(nfail<6){
		Bcalculated=false;
		CalcBTTrain();
		CalcSigmaA();
		CalcExactLogP();
		if(logP==logP && !isinf(logP) && logP>bestlogP){
			oldbestlogP=bestlogP;
			oldbestLambda=bestLambda;
			bestlogP=logP;
			bestLambda=LAMBDA;
			LAMBDA=bestLambda+dLambda;
		}
		else{
			bestLambda=oldbestLambda;
			bestlogP=oldbestlogP;
			nfail+=1;
			dLambda=dLambda/2.0;
			LAMBDA=bestLambda+dLambda;
			if(logP!=logP || isinf(logP)){
				CLog::Info("Note: LAMBDA set to "+to_string(LAMBDA)+" -- optimum value might be higher\n but determinant(B) cannot be calculated due to numerical accuracy problems.\n");
			}
		}
	}
	LAMBDA=bestLambda;
	CalcBTTrain();
	CalcSigmaA();
	CalcExactLogP();
}
