#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

void CSimplexSampler::GetC0DDprime(double LAMBDA,vector<double> &theta1,vector<double> &theta2,double &C0,double &D,double &Dprime){
	unsigned int ipar;
	double delThetaSquared,delTheta;
	delThetaSquared=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		delTheta=theta1[ipar]-theta2[ipar];
		delThetaSquared+=delTheta*delTheta;
	}
	C0=exp(-0.5*delThetaSquared/(LAMBDA*LAMBDA));
	D=C0*delThetaSquared;
	Dprime=D*delThetaSquared;
}

void CSimplexSampler::CalcIJK(double LAMBDA,double beta){
	unsigned int a,b,ipar;
	double lambda,gamma,dT2,tatb,ta2,tb2,Jab,Jba,Iab,Kab;
	double X,dXdgamma_a,dXdgamma_b,d2Xdgamma_adgamma_b;
	gamma=1.0/(LAMBDA*LAMBDA);
	lambda=2.0*gamma+beta;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			dT2=tatb=ta2=tb2=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				ta2+=pow(ThetaTrain[a][ipar],2);
				tb2+=pow(ThetaTrain[b][ipar],2);
				dT2+=pow(ThetaTrain[a][ipar]-ThetaTrain[b][ipar],2);
				tatb+=ThetaTrain[a][ipar]*ThetaTrain[b][ipar];
			}
			X=gamma*gamma*dT2+beta*gamma*(ta2+tb2);
			Iab=sqrt(pow(beta/lambda,NPars))*exp(-0.5*X/lambda);
			dXdgamma_a=gamma*dT2+beta*ta2;
			dXdgamma_b=gamma*dT2+beta*tb2;
			Jab=-(0.5*double(NPars)/lambda)*Iab+(0.5*X/(lambda*lambda))*Iab-(0.5*dXdgamma_a/lambda)*Iab;
			Jba=-(0.5*double(NPars)/lambda)*Iab+(0.5*X/(lambda*lambda))*Iab-(0.5*dXdgamma_b/lambda)*Iab;
			d2Xdgamma_adgamma_b=dT2;
			Kab=Jab*Jba/Iab +(0.5/(lambda*lambda))*Iab -(X/(lambda*lambda*lambda))*Iab
			+(0.5*(dXdgamma_a+dXdgamma_b)/(lambda*lambda))*Iab -(0.5*d2Xdgamma_adgamma_b/lambda)*Iab;
			Jab*=2.0;
			Jba*=2.0;
			Kab*=4.0;
			I(a,b)=Iab;
			J(a,b)=Jab+Jba;
			K(a,b)=Kab;
		}
	}
}

double CSimplexSampler::GetSigma2Bar(double LAMBDA,double ALPHA,double &detB){
	unsigned int a,b;
	Eigen::MatrixXd B,Binv,D,Dprime,B0,B2;
	double beta=3.0,SigmaA=1.0; // SigmaA shouldn't matter
	double dEdLambda2,d2logdetBdLambda2,Sigma2Bar;
	double bb,dd,ddprime,detB0,detB2,dLAMBDA=0.05;
	I.resize(NTrainingPts,NTrainingPts);
	J.resize(NTrainingPts,NTrainingPts);
	K.resize(NTrainingPts,NTrainingPts);
	CalcIJK(LAMBDA,beta);

	B.resize(NTrainingPts,NTrainingPts);
	Binv.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	Dprime.resize(NTrainingPts,NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			GetC0DDprime(LAMBDA,ThetaTrain[a],ThetaTrain[b],bb,dd,ddprime);
			B(a,b)=bb; D(a,b)=dd; Dprime(a,b)=ddprime;
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	Binv=B.inverse();

	//printf("|B|=%g\n",B.determinant());
	double Iterm,Jterm,Kterm;
	Iterm=(I*Binv*D*Binv*D*Binv).trace();
	Jterm=(J*Binv*D*Binv).trace();
	Kterm=(K*Binv).trace();
	dEdLambda2=Iterm+Jterm+Kterm;
	dEdLambda2=dEdLambda2/pow(LAMBDA,6);
	double S20=1.0-(I*Binv).trace();

	//dEdLambda2=(K*Binv).trace()+(J*Binv*D*Binv).trace()+(I*Binv*D*Binv*D*Binv).trace();
	if(dEdLambda2<-0.002){
		CLog::Info("dEdLambda2<0, ="+to_string(dEdLambda2)+"\n");
	}
	
	// Calculate d2|B|/dLambda^2
	detB=B.determinant();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			GetC0DDprime(LAMBDA-dLAMBDA,ThetaTrain[a],ThetaTrain[b],bb,dd,ddprime);
			B(a,b)=bb; D(a,b)=dd; Dprime(a,b)=ddprime;
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	detB0=B.determinant();
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			GetC0DDprime(LAMBDA+dLAMBDA,ThetaTrain[a],ThetaTrain[b],bb,dd,ddprime);
			B(a,b)=bb; D(a,b)=dd; Dprime(a,b)=ddprime;
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	detB2=B.determinant();

	d2logdetBdLambda2=(log(detB2)-2.0*log(detB)+log(detB0))/(dLAMBDA*dLAMBDA);
	if(d2logdetBdLambda2<0.0){
		CLog::Info("d2logdetBdLambda2="+to_string(d2logdetBdLambda2)+", is <0\n");
	}


	Eigen::Matrix2d W,Winv;

//-logZ=0.5*log(fabs(detB))+NTrainingPts*log(SigmaA)+0.5*y*Binv*y;

	W.resize(2,2);
	Winv.resize(2,2);
	Winv(0,0)=double(2*NTrainingPts)/(SigmaA*SigmaA);

	Winv(1,0)=Winv(0,1)=(-1.0/(SigmaA*pow(LAMBDA,3)))*((Binv*D).trace());

	Winv(1,1)=0.5*d2logdetBdLambda2 +(3.0/(2.0*pow(LAMBDA,4)))*(D*Binv).trace()
	+(1.0/pow(LAMBDA,6))*(D*Binv*D*Binv).trace() -(1.0/(2.0*pow(LAMBDA,6)))*(Dprime*Binv).trace();

	W=Winv.inverse();
	if(W(1,1)<0.0){
		printf("detB012=(%g,%g,%g)\n",detB0,detB,detB2);
		CLog::Info("W(1,1)<0, ="+to_string(W(1,1))+"\n");
		Sigma2Bar=1.0E99;
	}
	double S2_duetoLambda=dEdLambda2*W(1,1);
	Sigma2Bar=S20+S2_duetoLambda;
	//printf("Lambda=%7.4f, R=%7.4f: S20=%8.3f  S2_Lambda=%8.3f <Sigma_E^2>=%8.3f Iterm,Jterm,Kterm=(%8.3f,%8.3f,%8.3f), %g\n",LAMBDA,RGauss,S20,S2_duetoLambda=dEdLambda2*W(1,1),Sigma2Bar,Iterm,Jterm,Kterm,Iterm+Jterm+Kterm);
	return Sigma2Bar;

}