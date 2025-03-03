#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>
#include <algorithm>
#include "msu_smoothutils/randy.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

void CSimplexSampler::Optimize(double LAMBDASet,double ALPHAset){
	LAMBDA=LAMBDASet;
	ALPHA=ALPHAset;
	PLUS1=false;
	if(NMC==0)
		NMC=parmap.getI("Simplex_NMC",0);
	if(NMC==0){
		CLog::Info("Enter NMC: ");
		scanf("%u",&NMC);
	}
	if(OptimizeMethod=="MC"){
		PLUS1=false;
		if(NTrainingPts==0)
			NTrainingPts=parmap.getI("Simplex_NTrainingPts",0);
		if(NTrainingPts==0){
			CLog::Info("Enter NTrainingPts: ");
			scanf("%u",&NTrainingPts);
		}
		Optimize_MC();
	}
	else if(OptimizeMethod=="MCSphere"){
		PLUS1=true;
		NTrainingPts=parmap.getI("Simplex_NTrainingPts",0);
		if(NTrainingPts==0){
			CLog::Info("Enter NTrainingPts: ");
			scanf("%u",&NTrainingPts);
		}
		OptimizeSphere_MC();
	}
	else if(OptimizeMethod=="MCSimplex"){
		PLUS1=false;
		NTrainingPts=NPars+1;
		OptimizeSimplex_MC();
	}
	else if(OptimizeMethod=="MCSimplexPlus1"){
		PLUS1=true;
		NTrainingPts=NPars+2;
		OptimizeSimplex_MC();
	}
	else if(OptimizeMethod=="MCQuadratic"){
		PLUS1=true;
		NTrainingPts=(NPars+1)*(NPars+2)/2;
		Optimize_MC();
	}
	else if(OptimizeMethod=="MCSphereQuadratic"){
		PLUS1=true;
		NTrainingPts=(NPars+1)*(NPars+2)/2;
		OptimizeSphere_MC();
	}
}

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

double my_erfinv (double a)
{
	double p, r, t;
	t = fmaf (a, 0.0f - a, 1.0f);
	t = log(t);
	if (fabs(t) > 6.125f) { // maximum ulp error = 2.35793
		p =              3.03697567e-10; //  0x1.4deb44p-32 
		p = fma (p, t,  2.93243101e-8); //  0x1.f7c9aep-26 
		p = fma (p, t,  1.22150334e-6); //  0x1.47e512p-20 
		p = fma (p, t,  2.84108955e-5); //  0x1.dca7dep-16 
		p = fma (p, t,  3.93552968e-4); //  0x1.9cab92p-12 
		p = fma (p, t,  3.02698812e-3); //  0x1.8cc0dep-9 
		p = fma (p, t,  4.83185798e-3); //  0x1.3ca920p-8 
		p = fma (p, t, -2.64646143e-1); // -0x1.0eff66p-2 
		p = fma (p, t,  8.40016484e-1); //  0x1.ae16a4p-1 
	} else { // maximum ulp error = 2.35002
		p =              5.43877832e-9;  //  0x1.75c000p-28 
		p = fma (p, t,  1.43285448e-7); //  0x1.33b402p-23 
		p = fma (p, t,  1.22774793e-6); //  0x1.499232p-20 
		p = fma (p, t,  1.12963626e-7); //  0x1.e52cd2p-24 
		p = fma (p, t, -5.61530760e-5); // -0x1.d70bd0p-15 
		p = fma (p, t, -1.47697632e-4); // -0x1.35be90p-13 
		p = fma (p, t,  2.31468678e-3); //  0x1.2f6400p-9 
		p = fma (p, t,  1.15392581e-2); //  0x1.7a1e50p-7 
		p = fma (p, t, -2.32015476e-1); // -0x1.db2aeep-3 
		p = fma (p, t,  8.86226892e-1); //  0x1.c5bf88p-1 
	}
	r = a * p;
	return r;
}

void CSimplexSampler::SetThetaLatinHyperCube(vector<vector<double>> &theta){
	Crandy randy(time(NULL));
	double thetamax,dtheta,root2=sqrt(2.0);
	unsigned int ipar,itrain,is;
	theta.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++)
		theta[itrain].resize(NPars);
	
	vector<int> ishuffle(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ishuffle[itrain]=itrain;
	}
	
	for(ipar=0;ipar<NPars;ipar++){
		std::shuffle(std::begin(ishuffle), std::end(ishuffle), randy.mt);
		thetamax=priorinfo->ThetaPrior[ipar];
		dtheta=2.0*thetamax/double(NTrainingPts);
		for(itrain=0;itrain<NTrainingPts;itrain++){
			is=ishuffle[itrain];
			if(priorinfo->type[ipar]=="uniform"){
				theta[itrain][ipar]=-thetamax+dtheta*(is+randy.ran());
			}
			else if(priorinfo->type[ipar]=="gaussian"){
				dtheta=2.0/double(NTrainingPts);
				theta[itrain][ipar]=-1.0+dtheta*(is+randy.ran());
				theta[itrain][ipar]=my_erfinv(theta[itrain][ipar]);
				theta[itrain][ipar]*=thetamax*root2;
			}
		}
	}
	
	
}

void CSimplexSampler::Optimize_MC(){
	double Sigma2Bar,bestSigma2,dtheta,detB,W11,r,successrate;
	unsigned int imc,itrain,ipar,nfail=0,nsuccess=0;
	FILE *fptr,*fptr_vsNMC;
	Crandy randy(time(NULL));
	vector<vector<double>> besttheta;
	printf("NTrainingpts=%u, NMC=%u, LAMBDA=%g, ALPHA=%g\n",NTrainingPts,NMC,LAMBDA,ALPHA);
	dtheta=0.05/sqrt(double(NTrainingPts));
	
	/*
	double R0=1.0/sqrt(3.0);
	besttheta.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
	besttheta[itrain].resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
	besttheta[itrain][ipar]=R0*randy.ran_gauss();
	if(itrain==0)
	besttheta[itrain][ipar]=0.0;
	if(itrain==1 && ipar!=0){
	besttheta[itrain][ipar]=0.0;
	}
	}
	}
	*/

	
	if(FIRSTCALL==true){
		SetThetaLatinHyperCube(besttheta);
		SetThetaTrain(besttheta);
		FIRSTCALL=false;
	}
	else{
		printf("-------- check\n");
		besttheta=ThetaTrain;
	}
	
	Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
	bestSigma2=Sigma2Bar;
	printf("beginning: bestSigma2=%g\n",bestSigma2);
	
	string filename="Sigma2vsNMC/NTrain"+to_string(NTrainingPts)+".txt";
	fptr_vsNMC=fopen(filename.c_str(),"w");
	fprintf(fptr_vsNMC,"0 %g\n",bestSigma2);
	
	fptr=fopen("Sigma2vsNTrain_Latin.txt","a");
	fprintf(fptr,"%u %g\n",NTrainingPts,bestSigma2);
	fclose(fptr);
	
	for(imc=0;imc<NMC;imc++){
		for(itrain=0;itrain<NTrainingPts;itrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[itrain][ipar]=besttheta[itrain][ipar]+dtheta*randy.ran_gauss();
			}
		}
		Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
		//printf("Sigma2Bar=%g\n",Sigma2Bar);
		if(Sigma2Bar<bestSigma2){
			nsuccess+=1;
			bestSigma2=Sigma2Bar;
			for(itrain=0;itrain<NTrainingPts;itrain++){
				besttheta[itrain]=ThetaTrain[itrain];
			}
		}
		else
			nfail+=1;
		for(int n=0;n<100;n++){
			if(imc==pow(2,n)){
				fprintf(fptr_vsNMC,"%u %g\n",imc,bestSigma2);
				n+=1000;
			}
			if(imc<pow(2,n))
				n+=1000;
		}
		if((100*(imc+1)%NMC)==0){
			
			successrate=double(nsuccess)/double(nsuccess+nfail);
			printf("++++++++++ finished %g percent, bestSigma2=%g, success %%=%g, dtheta=%g ++++++++++\n",100.0*(imc+1.0)/double(NMC),bestSigma2,100.0*successrate,dtheta);
			dtheta*=0.05+4.0*successrate;
			nfail=nsuccess=0;
		}
		
	}
	printf("--- best Sigma2=%g ---\n",bestSigma2);
	fprintf(fptr_vsNMC,"%u %g\n",imc,bestSigma2);
	fclose(fptr_vsNMC);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		r=0.0;
		printf("%2d: ",itrain);
		for(ipar=0;ipar<NPars;ipar++){
			r+=pow(besttheta[itrain][ipar],2);
			printf("%8.5f ",besttheta[itrain][ipar]);
		}
		printf(": %8.5f\n",sqrt(r));
		
	}
	fptr=fopen("Sigma2vsNTrain.txt","a");
	fprintf(fptr,"%3d %8.6f\n",NTrainingPts,bestSigma2);
	fclose(fptr);
	
	fptr=fopen("Sigma2vsLambda.txt","a");
	fprintf(fptr,"%8.5f %8.6f\n",LAMBDA,bestSigma2);
	fclose(fptr);
	
	fptr=fopen("Sigma2vsNPars.txt","a");
	fprintf(fptr,"%8.5u %8.6f %u\n",NPars,bestSigma2,NTrainingPts);
	fclose(fptr);
	
	/*
	int jtrain;
	vector<vector<double>> ctheta;
	ctheta.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++)
	ctheta[itrain].resize(NTrainingPts);
	vector<double> rtrain(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
	rtrain[itrain]=0.0;
	for(ipar=0;ipar<NPars;ipar++){
	rtrain[itrain]+=besttheta[itrain][ipar]*besttheta[itrain][ipar];
	}
	rtrain[itrain]=sqrt(rtrain[itrain]);
	}
	for(ipar=0;ipar<NPars;ipar++)
	for(itrain=0;itrain<NTrainingPts;itrain++){
	for(jtrain=0;jtrain<NTrainingPts;jtrain++){
	double rij=0.0;
	if(rtrain[itrain]>0.1 && rtrain[jtrain]>0.1){
	ctheta[itrain][jtrain]=0.0;
	for(ipar=0;ipar<NPars;ipar++){
	rij+=pow(besttheta[itrain][ipar]-besttheta[jtrain][ipar],2);
	ctheta[itrain][jtrain]+=besttheta[itrain][ipar]*besttheta[jtrain][ipar];
	}
	rij=sqrt(rij);
	ctheta[itrain][jtrain]=ctheta[itrain][jtrain]/(rtrain[itrain]*rtrain[jtrain]);
	fprintf(fptr,"%8.5f ",ctheta[itrain][jtrain]);
	}
	}
	}			
	
	fptr=fopen("ctheta.txt","w");
	for(itrain=0;itrain<NTrainingPts;itrain++){
	sort(ctheta[itrain].begin(),ctheta[itrain].end());
	for(jtrain=0;jtrain<NTrainingPts;jtrain++){
	fprintf(fptr,"%8.5f ",ctheta[itrain][jtrain]);

	}
	fprintf(fptr,"\n");
	}
	fclose(fptr);
	*/
}

void CSimplexSampler::OptimizeSimplex_MC(){
	double R=1.2,bestSigma2,Sigma2bar,dR=0.1,W11,detB,bestR,successrate;
	unsigned int imc,nsuccess,nfail,itrain;
	vector<vector<double>> besttheta;
	Crandy randy(time(NULL));
	bestSigma2=1.0E99;
	bestR=R;
	dR=0.2/sqrt(double(NTrainingPts));
	besttheta.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		besttheta[itrain].resize(NPars);
	}
	if(PLUS1){
		SetThetaSimplexPlus1(R);
	}
	else
		SetThetaSimplex(R);

	Sigma2bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
	bestSigma2=Sigma2bar;
	printf("beginning: bestSigma2=%g\n",bestSigma2);
	
	nsuccess=nfail=0;
	for(imc=0;imc<NMC;imc++){
		R=bestR+dR*randy.ran_gauss();
		if(PLUS1)
			SetThetaSimplexPlus1(R);
		else
			SetThetaSimplex(R);
		Sigma2bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
		if(Sigma2bar<bestSigma2){
			bestSigma2=Sigma2bar;
			bestR=R;
			nsuccess+=1;
			for(itrain=0;itrain<NTrainingPts;itrain++)
				besttheta[itrain]=ThetaTrain[itrain];
		}
		else
			nfail+=1;
		if((100*(imc+1)%NMC)==0){
			successrate=double(nsuccess)/double(nsuccess+nfail);
			printf("+++++++++ finished %g percent, bestSigma2=%g, success %%=%g, dR=%g ++++++++\n",100.0*(imc+1.0)/double(NMC),bestSigma2,100.0*successrate,dR);
			dR*=0.05+4.0*successrate;
			nfail=nsuccess=0;
		}
	}

	printf("bestR=%g\n",bestR);
	SetThetaSimplex(bestR);
	SetThetaTrain(besttheta); 
	Sigma2bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
	bestSigma2=Sigma2bar;	
}

void CSimplexSampler::OptimizeSphere_MC(){
	double Sigma2Bar,bestSigma2=1.0E99,dtheta,detB,W11,r,R0=1.0,successrate;
	unsigned int imc,itrain,ipar,nfail=0,nsuccess=0;
	Crandy randy(time(NULL));
	vector<vector<double>> besttheta;

	NTrainingPts=(NPars+1)*(NPars+2)/2;
	dtheta=0.2/sqrt(double(NTrainingPts));
	besttheta.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		besttheta[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			besttheta[itrain][ipar]=R0*randy.ran_gauss();
			if(itrain==0)
				besttheta[itrain][ipar]=0.0;
			if(itrain==1 && ipar!=0){
				besttheta[itrain][ipar]=0.0;
			}
		}
	}
	for(itrain=1;itrain<NTrainingPts;itrain++){
		r=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			r+=besttheta[itrain][ipar]*besttheta[itrain][ipar];
		}
		r=sqrt(r);
		for(ipar=0;ipar<NPars;ipar++){
			besttheta[itrain][ipar]*=R0/r;
		}
	}
	SetThetaTrain(besttheta); 
	Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
	bestSigma2=Sigma2Bar;
	printf("beginning: bestSigma2=%g\n",bestSigma2);
	
	for(imc=0;imc<NMC;imc++){
		for(itrain=1;itrain<NTrainingPts;itrain++){
			r=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[itrain][ipar]=besttheta[itrain][ipar]+dtheta*randy.ran_gauss();
				if(itrain==1 && ipar!=0){
					ThetaTrain[itrain][ipar]=0.0;
				}
				r+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];	
			}
			r=sqrt(r);
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[itrain][ipar]*=ThetaTrain[1][0]/r;
			}
		}
		Sigma2Bar=GetSigma2Bar(LAMBDA,ALPHA,detB,W11);
		
		if(Sigma2Bar<bestSigma2){
			nsuccess+=1;
			bestSigma2=Sigma2Bar;
			//printf("%10d: bestSigma2=%g, R=%g, detB=%g, W11=%g\n",imc,bestSigma2,ThetaTrain[1][0],detB,W11);
			for(itrain=0;itrain<NTrainingPts;itrain++){
				besttheta[itrain]=ThetaTrain[itrain];
			}
		}
		else
			nfail+=1;
		if((100*(imc+1)%NMC)==0){
			successrate=double(nsuccess)/double(nsuccess+nfail);
			printf("++++++++ finished %g percent, bestSigma2=%g, bestR=%g, success %%=%g, dtheta=%g +++++++++++\n",100.0*(imc+1.0)/double(NMC),bestSigma2,fabs(besttheta[1][0]),100.0*successrate,dtheta);
			dtheta*=0.05+4.0*successrate;
			nfail=nsuccess=0;
		}
		
	}
	printf("--- best Sigma2=%g ---\n",bestSigma2);
	SetThetaTrain(besttheta);
}

double CSimplexSampler::GetSigma2Bar(double LAMBDA,double ALPHA,double &detB,double &W11){
	unsigned int a,b;
	Eigen::MatrixXd B,Binv,D,Dprime,BB0,BB1,BB2;
	double SigmaA=1.0; // SigmaA shouldn't matter
	double dEdLambda2,d2logdetBdLambda2,Sigma2Bar;
	double bb,dd,ddprime,detB0,detB1,detB2,dLAMBDA=0.01*LAMBDA,S20;
	double detfactor=4*sqrt(NTrainingPts);
	
	I.resize(NTrainingPts,NTrainingPts);
	J.resize(NTrainingPts,NTrainingPts);
	K.resize(NTrainingPts,NTrainingPts);
	CalcIJK(LAMBDA,priorinfo->ThetaPrior);
	
	/*
	printf("------ General  ------\n");
	for(a=0;a<NTrainingPts;a++){
	for(b=0;b<NTrainingPts;b++){
	printf("%10.3e ",K(a,b));
	}
	printf("\n");
	}
	CalcIJK_Gaussian(LAMBDA,1.0/sqrt(3.0));
	printf("------ Gaussian ------\n");
	for(a=0;a<NTrainingPts;a++){
	for(b=0;b<NTrainingPts;b++){
	printf("%10.3e ",K(a,b));
	}
	printf("\n");
	}
	*/

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
	S20=1.0-(I*Binv).trace();
	if(INCLUDE_LAMBDA_UNCERTAINTY){

		double Iterm,Jterm,Kterm;
		Iterm=(I*Binv*D*Binv*D*Binv).trace();
		Jterm=-(J*Binv*D*Binv).trace();
		Kterm=(K*Binv).trace();
		dEdLambda2=Iterm+Jterm+Kterm;
		dEdLambda2=dEdLambda2/pow(LAMBDA,6);

		if(dEdLambda2<-0.002){
			CLog::Info("dEdLambda2<0, ="+to_string(dEdLambda2)+"\n");
		}
	
		// Calculate d2|B|/dLambda^2
		BB1=B*detfactor;
		detB1=BB1.determinant();
	
		if(detB1!=detB1){
			CLog::Fatal("|B| != |B|\n");
		}
		if(fabs(detB1)<5.0E-324){
			CLog::Fatal("|B| too small, = "+to_string(detB1)+"\n");
		}
		for(a=0;a<NTrainingPts;a++){
			for(b=0;b<NTrainingPts;b++){
				GetC0DDprime(LAMBDA-dLAMBDA,ThetaTrain[a],ThetaTrain[b],bb,dd,ddprime);
				B(a,b)=bb; D(a,b)=dd; Dprime(a,b)=ddprime;
			}
		}
		for(a=0;a<NTrainingPts;a++)
			B(a,a)+=ALPHA*ALPHA;
		BB0=B*detfactor;
		detB0=BB0.determinant();
		for(a=0;a<NTrainingPts;a++){
			for(b=0;b<NTrainingPts;b++){
				GetC0DDprime(LAMBDA+dLAMBDA,ThetaTrain[a],ThetaTrain[b],bb,dd,ddprime);
				B(a,b)=bb; D(a,b)=dd; Dprime(a,b)=ddprime;
			}
		}
		for(a=0;a<NTrainingPts;a++)
			B(a,a)+=ALPHA*ALPHA;
		BB2=B*detfactor;
		detB2=BB2.determinant();

		d2logdetBdLambda2=(log(detB2)-2.0*log(detB1)+log(detB0))/(dLAMBDA*dLAMBDA);
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
		W11=W(1,1);
		if(W(1,1)<0.0){
			printf("detB012=(%g,%g,%g)\n",detB0,detB1,detB2);
			CLog::Info("W(1,1)<0, ="+to_string(W(1,1))+"\n");
			Sigma2Bar=1.0E99;
		}
		double S2_duetoLambda=dEdLambda2*W(1,1);
		//printf("S20=%g, S2_duetoLambda=%g\n",S20,S2_duetoLambda);
		Sigma2Bar=S20+S2_duetoLambda;
		detB=detB1;
		return Sigma2Bar;
	}
	else
		return S20;
}
