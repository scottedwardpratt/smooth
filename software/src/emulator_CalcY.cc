#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;

void CSmoothEmulator::CalcY(CModelParameters *modpars,double &Y,double &SigmaY_emulator){
	double y;
	Y=SigmaY_emulator=0.0;
	for(unsigned int isample=0;isample<NASample;isample++){
		y=smooth->CalcY(ASample[isample],LAMBDA,modpars->Theta);
		Y+=y;
		SigmaY_emulator+=y*y;
	}
	SigmaY_emulator=SigmaY_emulator/double(NASample);
	Y=Y/double(NASample);
	SigmaY_emulator=sqrt(fabs(SigmaY_emulator-Y*Y));
}

void CSmoothEmulator::CalcY(vector<double> Theta,double &Y,double &SigmaY_emulator){
	double y;
	Y=SigmaY_emulator=0.0;
	for(unsigned int isample=0;isample<NASample;isample++){
		y=smooth->CalcY(ASample[isample],LAMBDA,Theta);
		Y+=y;
		SigmaY_emulator+=y*y;
	}
	SigmaY_emulator=SigmaY_emulator/double(NASample);
	Y=Y/double(NASample);
	SigmaY_emulator=sqrt(fabs(SigmaY_emulator-Y*Y));
}
