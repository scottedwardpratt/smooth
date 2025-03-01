#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::GetYAndUncertainty(vector<double> &theta,double &Y,double &uncertainty){
	double unc2,unc2_Lambda; // squared uncertainty
	unsigned int a,b;
	vector<double> S;
	S.resize(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		S[a]=GetCorrelation(theta,smoothmaster->traininginfo->modelpars[a]->Theta);
		//printf("S[%u]=%g\n",a,S[a]);
	}

	Y=0.0;
	for(a=0;a<NTrainingPts;a++)
		Y+=chi[a]*S[a];


	unc2=GetCorrelation(theta,theta);
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2-=S[a]*Binv(a,b)*S[b];
		}
	}
	unc2*=SigmaA*SigmaA;
	
	//cout << Binv << endl;
	//printf("---- detB=%g, unc2=%g, SigmaA=%g\n",detB,unc2,SigmaA);
	
	if(INCLUDE_LAMBDA_UNCERTAINTY){
		unc2_Lambda=GetSigma2_Lambda(theta);
		//unc2_Lambda=0.0;
		unc2=unc2+unc2_Lambda;
	}
	//printf("unc2=%g, unc2_Lambda=%g, unc_Lambda/unc=%g\n",unc2,unc2_Lambda,sqrt(fabs(unc2_Lambda/unc2)));
	uncertainty=sqrt(fabs(unc2));
}
