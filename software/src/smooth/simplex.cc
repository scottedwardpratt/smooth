#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CSimplexSampler::CSimplexSampler(){
	randy=new Crandy(123);
	parmap.ReadParsFromFile("smooth_data/smooth_parameters/simplex_parameters.txt");
	string logfilename=parmap.getS("Simplex_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	TrainType=parmap.getI("Simplex_TrainType",1);
	string prior_info_filename="smooth_data/Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(prior_info_filename);
	CModelParameters::priorinfo=priorinfo;
	NPars=priorinfo->NModelPars;
	RGauss=parmap.getD("Simplex_RGauss",1.5);
	RGauss1=parmap.getD("Simplex_RGauss1",1.0);
	RGauss2=parmap.getD("Simplex_RGauss2",2.0);
}

void CSimplexSampler::SetThetaSimplex(){
	if(TrainType==1)
		SetThetaType1();
	else if(TrainType==2)
		SetThetaType2();
	else{
		CLog::Fatal("Inside CSimplexSampler::SetThetaSimplex, TrainType must be 1 or 2\n");
	}
	//CLog::Info("NTrainingPts="+to_string(NTrainingPts)+"\n");
}

void CSimplexSampler::SetThetaType1(){
	unsigned int ipar,itrain,jtrain;
	double R,z,RTrain;
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}

	RTrain=RGauss*sqrt(double(NPars)/3.0);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="gaussian")
				ThetaTrain[itrain][ipar]*=(RTrain/R);
			else{
				ThetaTrain[itrain][ipar]*=(RTrain/R);
			}
		}
	}

	double Rmax=0.95;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="uniform"){
				if(fabs(ThetaTrain[itrain][ipar])>Rmax){
					ThetaTrain[itrain][ipar]*=(Rmax/fabs(ThetaTrain[itrain][ipar]));
				}
			}
		}
	}
	
}

void CSimplexSampler::SetThetaType2(){
	unsigned int ipar,itrain,jtrain,N1,n;
	double R,z,RTrain1,RTrain2;
	
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	
	N1=NTrainingPts;
	n=N1;
	NTrainingPts+=N1*(N1-1)/2;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=N1;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
	}
	for(itrain=1;itrain<N1;itrain++){
		for(jtrain=0;jtrain<itrain;jtrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[n][ipar]=0.5*(ThetaTrain[itrain][ipar]+ThetaTrain[jtrain][ipar]);
			}
			n+=1;
		}
	}
	
	double R1,R2;
	itrain=0;
	R1=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		R1+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
	R1=sqrt(R1);
	itrain=NTrainingPts-1;
	R2=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		R2+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
	R2=sqrt(R2);
	RTrain1=RGauss1*sqrt(double(NPars)/3.0);
	RTrain2=RGauss2*sqrt(double(NPars)/3.0);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(itrain<=NPars){
				ThetaTrain[itrain][ipar]*=(RTrain1/R1);
			}
			else{
				ThetaTrain[itrain][ipar]*=(RTrain2/R2);
			}
		}
	}

	double Rmax=0.95;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="uniform"){
				if(fabs(ThetaTrain[itrain][ipar])>Rmax){
					ThetaTrain[itrain][ipar]*=(Rmax/fabs(ThetaTrain[itrain][ipar]));
				}
			}
		}
	}
	
}

void CSimplexSampler::WriteModelPars(){
	FILE *fptr;
	string filename,dirname,command;
	unsigned int itrain,ipar;
	vector<CModelParameters *> modelparameters(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelparameters[itrain]=new CModelParameters();
		for(ipar=0;ipar<NPars;ipar++){
			modelparameters[itrain]->Theta[ipar]=ThetaTrain[itrain][ipar];
		}
		modelparameters[itrain]->TranslateTheta_to_X();
	}
	for(itrain=0;itrain<NTrainingPts;itrain++){
		dirname="smooth_data/modelruns/run"+to_string(itrain);
		command="mkdir -p "+dirname;
		system(command.c_str());
		filename=dirname+"/mod_parameters.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%s %g\n",
			priorinfo->parname[ipar].c_str(),modelparameters[itrain]->X[ipar]);
		}
		fclose(fptr);
	}
}

void CSimplexSampler::GetSigmaBar(double LAMBDA,double ALPHA,double &SigmaBar2,double &detB,double &TrB,double &TrBinv){
	// assumes
	Eigen::MatrixXd B,Binv,Bminus,Bplus,D;
	Eigen::Matrix2d S,Sinv;
	vector<double> W;
	unsigned int a,b,ipar;
	double delTheta,delThetaSquared,Omega,exparg,beta,rbeta,thetabar;
	double dL=0.05,logdetB,logdetBplus,logdetBminus;
	W.resize(NPars);
	B.resize(NTrainingPts,NTrainingPts);
	Bminus.resize(NTrainingPts,NTrainingPts);
	Bplus.resize(NTrainingPts,NTrainingPts);
	Binv.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	beta=1.0+(2.0/(3.0*LAMBDA*LAMBDA));
	rbeta=1.0/sqrt(beta);

	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			delThetaSquared=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				delTheta=ThetaTrain[a][ipar]-ThetaTrain[b][ipar];
				delThetaSquared+=delTheta*delTheta;
			}
			B(a,b)=exp(-0.5*delThetaSquared/(LAMBDA*LAMBDA));
			D(a,b)=delThetaSquared*B(a,b);
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	Binv=B.inverse();

	double SigmaA=12.3456; // should be arbitrary (but check)
	S(0,0)=-2.0*NTrainingPts/(SigmaA*SigmaA);
	S(0,1)=S(1,0)=(Binv*D).trace()/SigmaA;
	S(1,1)=-S(0,1)*SigmaA;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			delThetaSquared=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				delTheta=ThetaTrain[a][ipar]-ThetaTrain[b][ipar];
				delThetaSquared+=delTheta*delTheta;
			}
			Bplus(a,b)=exp(-0.5*delThetaSquared/((LAMBDA+dL)*(LAMBDA+dL)));
			Bminus(a,b)=exp(-0.5*delThetaSquared/((LAMBDA-dL)*(LAMBDA-dL)));
		}
	}
	logdetB=log(B.determinant());
	logdetBminus=log(Bminus.determinant());
	logdetBplus=log(Bplus.determinant());
	S(1,1)-=0.5*(logdetBplus-2.0*logdetB+logdetBminus)/(dL*dL);
	Sinv=S.inverse();

	SigmaBar2=1.0;
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			Omega=1.0;
			for(ipar=0;ipar<NPars;ipar++){
				thetabar=0.5*(ThetaTrain[a][ipar]+ThetaTrain[b][ipar]);
				exparg=(-0.5/(LAMBDA*LAMBDA))
				*(ThetaTrain[a][ipar]*ThetaTrain[a][ipar]+ThetaTrain[b][ipar]*ThetaTrain[b][ipar])
				+thetabar*thetabar/(1.5*beta*pow(LAMBDA,4));
				double wiab=rbeta*exp(exparg);
				Omega*=wiab;
				if(rbeta*exp(exparg)>1.0){
					CLog::Fatal("exp(exparg)/sqrt(beta)="+to_string(exp(exparg)/sqrt(beta))+"\n");
				}
			}
			SigmaBar2-=Omega*Binv(a,b);
		}
	}

	SigmaBar2-+Sinv(1,1)*.....;



	detB=B.determinant();
	TrB=B.trace();
	TrBinv=Binv.trace();

}

void CSimplexSampler::GetC0(vector<double> &theta1,vector<double> &theta2,double &C0,double){
	double answer,delTheta,delThetaSquared;
	int a,b;
	
	delThetaSquared=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		delTheta=ThetaTrain[a][ipar]-ThetaTrain[b][ipar];
		delThetaSquared+=delTheta*delTheta;
	}
	return exp(-0.5*delThetaSquared/(LAMBDA*LAMBDA));

}

void CSimplexSampler::GetC0D0(vector<double> &theta1,vector<double> &theta2,double &C0,double &D0){
	int ipar;
	double delThetaSquared,delTheta;
	delThetaSquared=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		delTheta=ThetaTrain[a][ipar]-ThetaTrain[b][ipar];
		delThetaSquared+=delTheta*delTheta;
	}
	C0=exp(-0.5*delThetaSquared/(LAMBDA*LAMBDA));
	D0=C0*delThetaSquared;
}

double CSimplexSampler::GetSigma2(double LAMBDA,double ALPHA,vector<double> &theta){
	unsigned int ipar,a,b;
	Eigen::MatrixXd B,Binv,D,Dinv;
	double sigma2,L2;
	double dthetaa,dthetab,dthetasuma,dthetasumb;
	double alpha,beta,gamma;

	B.resize(NTrainingPts,NTrainingPts);
	Binv.resize(NTrainingPts,NTrainingPts);
	D.resize(NTrainingPts,NTrainingPts);
	Dinv.resize(NTrainingPts,NTrainingPts);
	for(a=0;i<NPars;a++){
		for(b=0;b<NPars;b++){
			GetC0D0(ThetaTrain[a],ThetaTrain[b],B[a][b],D[a][b]);
		}
	}
	for(a=0;a<NTrainingPts;a++)
		B(a,a)+=ALPHA*ALPHA;
	Binv=B.inverse();
	Dinv=D.inverse();

	alpha=1.0/sqrt(3.0);
	L2=LAMBDA*LAMBDA;
	beta=0.5*L2*alpha*alpha/(alpha*alpha+0.5*L2);


	Eigen::MatrixXD I0,I1a,I1b,I2ab;
	double F0;
	I0.resize(NTrainingPts,NTrainingPts);
	I1a.resize(NTrainingPts,NTrainingPts);
	I1b.resize(NTrainingPts,NTrainingPts);
	I2ab.resize(NTrainingPts,NTrainingPts);

	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			I0(a,b)=I1a(a,b)=I1b(a,b)=I2ab(a,b)=1.0;
			for(ipar=0;ipar<NPars;ipar++){
				gamma=(ThetaTrain[a]-ThetaTrain[b])*beta*beta/L2;
				F0=(beta/gamma)
				*exp(-0.5*(pow(ThetaTrain[a][ipar],2)+pow(ThetaTrain[b][ipar],2))/L2)
				*exp(0.5*pow(beta/LAMBDA,2)*pow(ThetaTrain[a][ipar]-ThetaTrain[b][ipar],2));
				I0(a,b)*=F0;
				I1a(a,b)*=F0*(beta*beta+pow(gamma-ThetaTrain[a][ipar],2));
				I1b(a,b)*=F0*(beta*beta+pow(gamma-ThetaTrain[b][ipar],2));
				I2ab(a,b)*=F0*(3.0*pow(beta,4)+4.0*(gamma-ThetaTrain[a][ipar])*(gamma-ThetaTrain[b][ipar])*beta*beta
					+pow(gamma-ThetaTrain[a][ipar],2)*pow(gamma-ThetaTrain[a][ipar],2));
			}
		}
	}
	double dEdLambda2=(I2ab*Binv).trace();
	dEdLambda2-=(I1a*Binv*D*Binv).trace();
	dEdLambda2-=(I1b*Binv*D*Binv).trace();
	dEdLambda2+=(I0*Binv*D*Binv*D*Binv).trace();

	Eigen::Matrix2d W,Winv;
	W(0,0)=-2.0*NTrainingPts/(SigmaA*SigmaA);
	W(1,0)=W(0,1)=((Binv*D)/trace())/SigmaA;
	W(1,1)=-0.5*(Binv*D).trace()+0.5*d2logdetBdLambda2;
	Winv=W.inverse();

	Sigma2=(I0*Binv).trace()+dEdLambda2/Winv(2,2);

}
