#include "emulator.h"
#include "smooth.h"
using namespace std;



void CSimplexSampler::SetThetaSimplex(vector<vector<double>> &ThetaTrain,unsigned int &NTrain){
	if(TrainRank==1)
		SetThetaRank1(ThetaTrain,NTrain);
	else if(TrainRank==2)
		SetThetaRank2(ThetaTrain,NTrain);
	else{
		printf("Smooth_ThetaRank in parmap must be 0 or 1\n");
		exit(1);
	}
}

void CSimplexSampler::SetThetaRank1(vector<vector<double>> &ThetaTrain,unsigned int &NTrain){
	unsigned int ipar,itrain,jtrain;
	double R,z,RTrain;
	RTrain=1.0-1.0/double(NPars+1);
	NTrain=NPars+1;
	ThetaTrain.resize(NTrain);
	for(itrain=0;itrain<NTrain;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrain;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	
	for(itrain=0;itrain<NTrain;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]*=(RTrain/R);
		}
	}
}

void CSimplexSampler::SetThetaRank2(vector<vector<double>> &ThetaTrain, unsigned int &NTrain){
	unsigned int ipar,itrain,jtrain,N1,n;
	double R,z,RTrain;
	RTrain=1.0-1.0/double(NPars+1);
	NTrain=NPars+1;
	ThetaTrain.resize(NTrain);
	for(itrain=0;itrain<NTrain;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrain;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	
	N1=NTrain;
	n=N1;
	NTrain+=N1*(N1-1)/2;
	ThetaTrain.resize(NTrain);
	for(itrain=N1;itrain<NTrain;itrain++){
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
	
	for(itrain=0;itrain<NTrain;itrain++){
		if(itrain==0){
			R=0.0;
			for(ipar=0;ipar<NPars;ipar++)
				R+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
			R=sqrt(R);
		}
		for(ipar=0;ipar<=NPars;ipar++){
			ThetaTrain[itrain][ipar]*=(RTrain/R);
		}
	}
}
