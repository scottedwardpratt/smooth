#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

CLLCalc::CLLCalc(CSmoothMaster *master_set){
	master=master_set;
	NPars=master->NPars;
	priorinfo=master->priorinfo;
	obsinfo=master->obsinfo;
	obsinfo->ReadExperimentalInfo("Info/experimental_info.txt"); // might want to change this later to be more flexible
	NObs=obsinfo->NObservables;
	Y.resize(NObs);
	dYdTheta.resize(NObs);
	SigmaY.resize(NObs);
	SigmaY_emulator.resize(NObs);
	for(unsigned int iy=0;iy<NObs;iy++){
		dYdTheta[iy].resize(NPars);
	}	
}

void CLLCalc::CalcLL(CModelParameters *modpars,double &LL){
	LL=-1.0;
}

void CLLCalc::CalcLLPlusDerivatives(CModelParInfo *modpars,double &LL,vector<double> &dLL_dtheta){
	LL=1000.0;
}

void CLLCalcSmooth::CalcLL(CModelParameters *modpars,double &LL){
	unsigned int iobs;
	master->CalcAllY(modpars,Y,SigmaY_emulator);
	LL=0.0;
	for(iobs=0;iobs<NObs;iobs++){
		sigma2=Sigmay_emulator[iY]*SigmaY_emulator[iy]+obsinfo->SigmaExp=[iy]*obsinfo->SigmaExp[iy];
		LL+=0.5*pow(Y[iobs]-obsinfo->Yexp[iobs],2)/sigma2;
	}
	
}

void CLLCalcSmooth::CalcLLPlusDerivatives(CModelParameters *modpars,double &LL){
	unsigned int iobs;
	master->CalcAllY(modpars,Y,SigmaY_emulator);
	LL=0.0;
	for(iobs=0;iobs<NObs;iobs++){
		sigma2=Sigmay_emulator[iY]*SigmaY_emulator[iy]+obsinfo->SigmaExp=[iy]*obsinfo->SigmaExp[iy];
		LL+=0.5*pow(Y[iobs]-obsinfo->Yexp[iobs],2)/sigma2;
	}
	
}