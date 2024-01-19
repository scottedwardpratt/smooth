#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;

void CSmoothEmulator::CalcY(CModelParameters *modpars,double &Y,double &SigmaY_emulator){
	CalcY(modpars->Theta,Y,SigmaY_emulator);
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


void CSmoothEmulator::CalcYDYDTheta(CModelParameters *modpars,double &Y,vector<double> &dYdTheta,double &SigmaY){
	CalcYDYDTheta(modpars->Theta,Y,dYdTheta,SigmaY);
}

void CSmoothEmulator::CalcYDYDTheta(vector<double> Theta,double &Y,vector<double> &dYdTheta,double &SigmaY){
	double y;
	unsigned int ipar;
	vector<double> dydtheta;
	dYdTheta.resize(NPars);
	dydtheta.resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
		dYdTheta[ipar]=0.0;
	}
	Y=SigmaY=0.0;
	for(unsigned int isample=0;isample<NASample;isample++){
		smooth->CalcYDYDTheta(ASample[isample],LAMBDA,Theta,y,dydtheta);
		Y+=y/double(NASample);
		SigmaY+=y*y/double(NASample);
		for(ipar=0;ipar<NPars;ipar++){
			dYdTheta[ipar]+=dydtheta[ipar]/double(NASample);
		}
	}
	SigmaY=sqrt(fabs(SigmaY-Y*Y));
}


