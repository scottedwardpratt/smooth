#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

double PI=4.0*atan(1.0);

/*
void CSimplexSampler::CalcIJK(double Lambda,double Alpha,vector<double> &Rprior){
	unsigned int a,b,ipar;
	vector<double> I,Jfact,Jafact,Jbfact,Kabfact;
	I.resize(NPars);
	Jfact.resize(NPars);
	Jafact.resize(NPars);
	Jbfact.resize(NPars);
	Kabfact.resize(NPars);
	
	
	
	
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			for(ipar=0;ipar<NPars;ipar++){
				if(priorinfo->type[ipar]=="uniform"){
					GetIiJiKiUniform(priorinfo->Rprior[ipar],Lambda,ThetaTrain[a][ipar],ThetaTrain[b][ipar],
					I[ipar],Jafact[ipar],Jbfact[ipar],Kabfact[ipar]);
				}
				else if(priorinfo->type[ipar]=="gaussian"){
					GetIiJiKiGaussian(priorinfo->Rprior[ipar],Lambda,ThetaTrain[a][ipar],ThetaTrain[b][ipar],
					I[ipar],Jafact[ipar],Jbfact[ipar],Kabfact[ipar]);
				}
			}
			

		else if(priorinfo->type[ipar]=="uniform"){
			CLog::Info("not ready\n");
		}
		else
			CLog::Fatal("priorinfo->type not gaussian or uniform\n");
	}
}
*/

void CSimplexSampler::CalcIJK_Gaussian(double LAMBDA,double beta){ // only works when all have gaussian priors with same beta
	unsigned int a,b,ipar;
	double lambda,gamma,dT2,tatb,ta2,tb2,Jab,Jba,Iab,Kab;
	double X,dXdgamma_a,dXdgamma_b,d2Xdgamma_adgamma_b;
	gamma=1.0/(LAMBDA*LAMBDA);
	lambda=2.0*gamma+3*beta;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			dT2=tatb=ta2=tb2=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				ta2+=pow(ThetaTrain[a][ipar],2);
				tb2+=pow(ThetaTrain[b][ipar],2);
				dT2+=pow(ThetaTrain[a][ipar]-ThetaTrain[b][ipar],2);
				tatb+=ThetaTrain[a][ipar]*ThetaTrain[b][ipar];
			}
			X=gamma*gamma*dT2+3*beta*gamma*(ta2+tb2);
			Iab=sqrt(pow(3*beta/lambda,NPars))*exp(-0.5*X/lambda);
			dXdgamma_a=gamma*dT2+3*beta*ta2;
			dXdgamma_b=gamma*dT2+3*beta*tb2;
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

void CSimplexSampler::GetIiJiKiGaussian(double Rprior,double Lambda,double theta_a,double theta_b,
double &I,double &Jaterm,double &Jbterm,double &Kabterm){
	double X,gamma,alpha,lambda,deltheta2,sumt2,Jterm;
	gamma=1.0/(Lambda*Lambda);
	alpha=1.0/(Rprior*Rprior);
	lambda=2.0*gamma+alpha;
	X=gamma*gamma*pow(theta_a-theta_b,2)+alpha*gamma*(theta_a*theta_a+theta_b*theta_b);
	deltheta2=(theta_a-theta_b)*(theta_a-theta_b);
	sumt2=theta_a*theta_a+theta_b*theta_b;
	
	I=sqrt(alpha/lambda)*exp(-0.5*X/lambda);
	Jterm=(-1.0/lambda)
		-deltheta2*((gamma*gamma+alpha*gamma)/(lambda*lambda))
		+sumt2*(-0.5*alpha*alpha/(lambda*lambda));	
	Jaterm=0.5*Jterm-0.25*(alpha/lambda)*(theta_a*theta_a-theta_b*theta_b);
	Jbterm=Jterm-Jaterm;
	//Jab=I*Jterm;
	Kabterm=Jaterm*Jbterm+0.5/(lambda*lambda)
		+(0.5*alpha*alpha*sumt2-0.5*(2*gamma*gamma+alpha*lambda)*deltheta2)/pow(lambda,3);
	Kabterm*=4;
	//Kab=4*I*Kterm;
	//Jab*=-2;
	Jaterm*=-2;
	Jbterm*=-2;
	
}

void CSimplexSampler::GetIiJiKiUniform(double Rprior,double Lambda,double theta_a,double theta_b,
double &I,double &Jaterm,double &Jbterm,double &Kabterm){
	double Ja,Jb,J,Kab;
	double deltheta,thetabar,rootgamma,Xplus,Xminus,P,Y,W,bplus2,bminus2,deltheta2,bplus,bminus;
	double gamma=1.0/(Lambda*Lambda);
	rootgamma=sqrt(gamma);
	thetabar=0.5*(theta_a+theta_b);
	deltheta=0.5*(theta_a-theta_b);
	deltheta2=deltheta*deltheta;
	bplus=Rprior-thetabar;
	bminus=-Rprior-thetabar;
	bplus2=bplus*bplus;
	bminus2=bminus*bminus;
	Xplus=exp(-gamma*bplus2);
	Xminus=exp(-gamma*bminus2);
	P=(0.5/Rprior)*exp(-gamma*deltheta*deltheta);
	W=(1.0/rootgamma)*(sqrt(PI)/2.0)*(erf(rootgamma*(Rprior-thetabar))-erf(rootgamma*(-Rprior-thetabar)));
	I=P*W;
	
	J=-(0.5/gamma)*I-deltheta2*I;
	Y=bplus*Xplus-bminus*Xminus;
	J+=(0.5/gamma)*P*Y;
	Ja=J+P*Xplus*deltheta/gamma-P*Xminus*deltheta/gamma;
	Jb=2.0*J-Ja;	

	Kab=J*((-0.5/gamma)-deltheta2);
	Kab+=0.5*I/(gamma*gamma);
	Kab=Kab-P*Xplus*bplus*((0.5/(gamma*gamma))+0.5*(deltheta2/gamma)+0.5*bplus*bplus/gamma);
	Kab=Kab+P*Xminus*bminus*((0.5/(gamma*gamma))+0.5*(deltheta2/gamma)+0.5*bminus*bminus/gamma);
	
	Kab-=I*(2.0*deltheta2/gamma);
	Kab=Kab+P*Xplus*bplus*2.0*deltheta2/gamma;
	Kab=Kab-P*Xminus*bminus*2.0*deltheta2/gamma;

	J=J*2.0;
	Jaterm=Ja/I;
	Jbterm=Jb/I;
	Kabterm=Kab/I;

}