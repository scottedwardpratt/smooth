#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUUtils;

void CSmoothEmulator::GetYAndUncertainty(vector<double> &Theta_s,double &Y,double &uncertainty){
	double unc2; // squared uncertainty
	unsigned int a,b;
	vector<double> S;
	printf("ALPHA=%g\n",ALPHA);
	S.resize(NTrainingPts);
	for(a=0;a<NTrainingPts;a++){
		S[a]=GetCorrelation(Theta_s,smoothmaster->traininginfo->modelpars[a]->Theta);
	}

	Y=0.0;
	for(a=0;a<NTrainingPts;a++)
		Y+=chi[a]*S[a];

	unc2=GetCorrelation(Theta_s,Theta_s);
	for(a=0;a<NTrainingPts;a++){
		for(b=0;b<NTrainingPts;b++){
			unc2-=S[a]*Binv(a,b)*S[b];
		}
	}
	uncertainty=SigmaA*sqrt(fabs(unc2));
}
