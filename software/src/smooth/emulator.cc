#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

unsigned int CSmoothEmulator::NPars=0;
CSmooth *CSmoothEmulator::smooth=NULL;
CSmoothMaster *CSmoothEmulator::smoothmaster=NULL;
CparameterMap *CSmoothEmulator::parmap=NULL;
Crandy *CSmoothEmulator::randy=NULL;
unsigned int CSmoothEmulator::NTrainingPts=0;

CSmoothEmulator::CSmoothEmulator(string observable_name_set,bool pca_ignore_set){
	observable_name=observable_name_set;
	NTrainingPts=smoothmaster->traininginfo->NTrainingPts;
	LAMBDA=parmap->getD("SmoothEmulator_LAMBDA",2.0);
	ConstrainA0=parmap->getB("SmoothEmulator_ConstrainA0",false);
	iY=smoothmaster->observableinfo->GetIPosition(observable_name);
	pca_ignore=pca_ignore_set;
	ThetaTrain.clear();
	TTrain.clear();
	Bcalculated=false;
	FixSigmaA=false;
}

void CSmoothEmulator::CalcBTTrain(){
	if(!Bcalculated){
		unsigned int a,b,ic,NCoefficients=smooth->NCoefficients;
		ThetaTrain.clear();
		TTrain.clear();
		ThetaTrain.resize(NTrainingPts);
		TTrain.resize(NTrainingPts);
		for(unsigned int itrain=0;itrain<NTrainingPts;itrain++){
			ThetaTrain[itrain].resize(NPars);
			TTrain[itrain].resize(NCoefficients);
			ThetaTrain[itrain]=smoothmaster->traininginfo->modelpars[itrain]->Theta;
			for(ic=0;ic<smooth->NCoefficients;ic++){
				TTrain[itrain][ic]=smooth->GetT(ic,LAMBDA,ThetaTrain[itrain]);
			}
		}
		B.resize(NTrainingPts,NTrainingPts);
		Binv.resize(NTrainingPts,NTrainingPts);
		for(a=0;a<NTrainingPts;a++){
			for(b=0;b<NTrainingPts;b++){
				B(a,b)=0.0;
				for(ic=0;ic<NCoefficients;ic++){
					B(a,b)+=TTrain[a][ic]*TTrain[b][ic];
				}
			}
		}
		B.resize(NTrainingPts,NTrainingPts);
		Binv.resize(NTrainingPts,NTrainingPts);
		for(a=0;a<NTrainingPts;a++){
			for(b=0;b<NTrainingPts;b++){
				B(a,b)=0.0;
				for(ic=0;ic<NCoefficients;ic++){
					B(a,b)+=TTrain[a][ic]*TTrain[b][ic];
				}
			}
		}
		
		if(!FixSigmaA)
			CalcSigmaA();  // note this ignores effect of model randomness(SigmaYTrain) in setting SigmaA
		//printf("---- SigmaA=%g\n",SigmaA);
	
		for(a=0;a<NTrainingPts;a++){
			B(a,a)+=pow(smoothmaster->traininginfo->SigmaYTrain[iY][a]/SigmaA,2);
		}
		Binv=B.inverse();
		Bcalculated=true;
	}
}

void CSmoothEmulator::Tune(){
	unsigned int a,b,ic,NCoefficients=smooth->NCoefficients;
	vector<double> chi;
	
	ABest.resize(NCoefficients);
	CalcBTTrain();
	
	chi.resize(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		chi[a]=0.0;
		for(b=0;b<NTrainingPts;b++){
			chi[a]+=Binv(a,b)*smoothmaster->traininginfo->YTrain[iY][b];
		}
	}
	for(ic=0;ic<NCoefficients;ic++){
		ABest[ic]=0.0;
		for(a=0;a<NTrainingPts;a++){
			ABest[ic]+=chi[a]*TTrain[a][ic];
		}
	}
	//CalcExactLogP();
}

void CSmoothEmulator::TuneForSigmaA(double SigmaAset){
	SetSigmaA(SigmaAset);
	Tune();
}

void CSmoothEmulator::SetA_Zero(vector<double> &A){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=0.0;
	}
}

void CSmoothEmulator::SetA_RanGauss(double SigmaA,vector<double> &A){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=SigmaA*randy->ran_gauss();
	}
}

void CSmoothEmulator::SetA_Constant(double SigmaA,vector<double> &A){
	if(!pca_ignore){
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=SigmaA;
	}
}

void CSmoothEmulator::SetA_RanSech(double SigmaA,vector<double> &A){
	if(!pca_ignore){
		const double PI=4.0*atan(1.0);
		double r=1.0-2.0*randy->ran();
		for(unsigned int ic=0;ic<A.size();ic++)
			A[ic]=SigmaA*2.0*atanh(tan(0.5*PI*(r-0.5)));
	}
}
