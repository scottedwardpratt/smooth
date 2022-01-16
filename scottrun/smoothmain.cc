#include "parametermap.h"
#include "constants.h"
#include "smooth.h"
#include "emulator.h"
#include "gslmatrix.h"

using namespace std;
int main(int argc,char *argv[]){
	double LAMBDA=10.0,AMAG=1.5;
	vector<double> realA,A;
	vector<vector<double>> ASample;
	unsigned int nsample=10,NMC=10000,isample,ipar,itest;

	unsigned int NPars,NTrainingPts,iTrain;
	printf("Enter NPars: ");
	scanf("%u",&NPars);
	vector<double> Theta,Lambda;
	vector<double> YTrain;
	vector<vector<double>> ThetaTrain;
	Theta.resize(NPars);
	CEmulator emulator(NPars);
	emulator.SetSize(A,Lambda);
	emulator.SetLambda_Constant(LAMBDA,Lambda);
	emulator.SetThetaRank1(NTrainingPts,ThetaTrain);
	printf("NTrainingPts=%u\n",NTrainingPts);

	YTrain.resize(NTrainingPts);
	realA=A;

	emulator.SetA_RanGauss(AMAG,realA);
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		YTrain[iTrain]=emulator.smooth->CalcY(realA,Lambda,ThetaTrain[iTrain]);
	}

	//emulator.CalcAFromTraining(ThetaTrain,YTrain,A,Lambda);
	emulator.TuneA(ThetaTrain,YTrain,AMAG,NMC,A,Lambda);
	printf("Test at training points\n");
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		printf("%g =? %g\n",emulator.smooth->CalcY(A,Lambda,ThetaTrain[iTrain]),
			YTrain[iTrain]);
	}

	ASample.resize(nsample);
	double yreal,y,accuracy=0.0;
	unsigned int ntest=10;
	for(isample=0;isample<nsample;isample++){
		//emulator.SetA_RanGauss(AMAG,realA);
		//for(iTrain=0;iTrain<NTrainingPts;iTrain++){
			//YTrain[iTrain]=emulator.smooth->CalcY(realA,Lambda,ThetaTrain[iTrain]);
		//}
		emulator.TuneA(ThetaTrain,YTrain,AMAG,NMC,A,Lambda);
		/*
		printf("Test at training points\n");
		for(iTrain=0;iTrain<NTrainingPts;iTrain++){
			printf("%g =? %g\n",emulator.smooth->CalcY(A,Lambda,ThetaTrain[iTrain]),
				YTrain[iTrain]);
		}
		*/
		printf("Now test away from training points\n");
		for(itest=0;itest<ntest;itest++){
			for(ipar=0;ipar<NPars;ipar++){
				Theta[ipar]=1.0-2.0*emulator.randy->ran();
			}
			y=emulator.smooth->CalcY(A,Lambda,Theta);
			yreal=emulator.smooth->CalcY(realA,Lambda,Theta);
			accuracy+=(y-yreal)*(y-yreal);
			//printf("%g =? %g\n",y,yreal);
		}
		ASample[isample]=A;
	}
	accuracy=accuracy/double(nsample*ntest);
	printf("accuracy=%g\n",sqrt(accuracy));
	
	return 0;
}