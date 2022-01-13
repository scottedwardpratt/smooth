#include "parametermap.h"
#include "constants.h"
#include "smooth.h"

using namespace std;
int main(int argc,char *argv[]){
	double F,Fbar=0.0;
	double LAMBDA=5.0;
	vector<double> theta;
	unsigned int maxrank=5,NPars=20,ntries=1,itry,ipar;
	CRandy *randy=new CRandy(-time(NULL));
	theta.resize(NPars);
	
	CSmooth smooth(maxrank,NPars);
	smooth.randy=randy;
	smooth.SetLambda(LAMBDA);
	smooth.SetA_Constant(1.0);
	
	Fbar=0.0;
	for(itry=0;itry<ntries;itry++){
		for(ipar=0;ipar<NPars;ipar++){
			theta[ipar]=1.0;
			//theta[ipar]=randy->ran_gauss();
		}
		F=smooth.CalcY(theta);
		Fbar+=F;
		//printf("F=%g\n",F);
	}
	Fbar=F/double(ntries);
	printf("<F>=%g\n",F);
	

	unsigned int iTrain;
	unsigned int NTrainingPts;
	printf("How many TrainingPts?    NPars=%u, For all quadratic, NTrain=%u\n",NPars,NPars+1+NPars*(NPars+1)/2);
	scanf("%u",&NTrainingPts);
	
	CSmooth realsmooth(maxrank,NPars);
	realsmooth.randy=randy;
	realsmooth.SetLambda(LAMBDA);
	realsmooth.SetA_RanGauss(1.0);
	
	vector<vector<double>> thetatheta;
	vector<double> YTrain;
	YTrain.resize(NTrainingPts);
	thetatheta.resize(NTrainingPts);
	for(iTrain=0;iTrain<NTrainingPts;iTrain++){
		thetatheta[iTrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			thetatheta[iTrain][ipar]=1.0-2.0*randy->ran();
		}
		YTrain[iTrain]=realsmooth.CalcY(thetatheta[iTrain]);
	}
	
	double Yfit,Yreal;
	CSmooth smooth_fit;
	smooth_fit.Copy(&smooth);
	for(itry=0;itry<10;itry++){
		printf("---------------------------------\n");
		printf("-- Check traininig pts\n");
		smooth_fit.SetA_RanGauss(0.1);
		smooth_fit.CalcAFromTraining(thetatheta,YTrain);
		for(iTrain=0;iTrain<NTrainingPts;iTrain++){
			Yfit=smooth_fit.CalcY(thetatheta[iTrain]);
			printf("%g =? %g\n",Yfit,YTrain[iTrain]);
		}
		printf("-- Test predictions\n");
		for(unsigned int itheta=0;itheta<10;itheta++){
			for(ipar=0;ipar<NPars;ipar++)
				theta[ipar]=1.0-2.0*randy->ran();
			Yfit=smooth_fit.CalcY(theta);
			Yreal=realsmooth.CalcY(theta);
			printf("%g =? %g\n",Yfit,Yreal);
		}
		
		
	}
	return 0;
}