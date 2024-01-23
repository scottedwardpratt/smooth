#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

CLLCalc::CLLCalc(){
	bestLL=-1.0E100;
}

CLLCalc::CLLCalc(CSmoothMaster *master_set){
	master=master_set;
	NPars=master->NPars;
	priorinfo=master->priorinfo;
	obsinfo=master->observableinfo;
	obsinfo->ReadExperimentalInfo("Info/experimental_info.txt"); // might want to change this later to be more flexible
	NObs=obsinfo->NObservables;
	Y.resize(NObs);
	dYdTheta.resize(NObs);
	SigmaY.resize(NObs);
	SigmaY_emulator.resize(NObs);
	for(unsigned int iy=0;iy<NObs;iy++){
		dYdTheta[iy].resize(NPars);
	}	
	bestLL=-1.0E100;
}

CLLCalcSmooth::CLLCalcSmooth(CSmoothMaster *master_set){
	master=master_set;
	NPars=master->NPars;
	priorinfo=master->priorinfo;
	obsinfo=master->observableinfo;
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
	(void) modpars;
	(void) LL;
}

void CLLCalc::CalcLLPlusDerivatives(CModelParameters *modpars,double &LL,vector<double> &dLL_dtheta){
	(void) modpars;
	(void) LL;
	(void) dLL_dtheta;
}

void CLLCalcSmooth::CalcLL(CModelParameters *modpars,double &LL){
	unsigned int iy,ipar;
	double sigma2;
	bool insidebounds=true;
	LL=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		if(modpars->priorinfo->type[ipar]=="uniform" && fabs(modpars->Theta[ipar])>1.0){
			LL=-1.0E100;
			insidebounds=false;
		}
	}
	if(insidebounds){
		master->CalcAllY(modpars,Y,SigmaY_emulator);
		for(iy=0;iy<NObs;iy++){
			sigma2=SigmaY_emulator[iy]*SigmaY_emulator[iy]+obsinfo->SigmaExp[iy]*obsinfo->SigmaExp[iy];
			LL-=0.5*pow(Y[iy]-obsinfo->YExp[iy],2)/sigma2;
		}
		for(ipar=0;ipar<NPars;ipar++){
			if(modpars->priorinfo->type[ipar]=="gaussian"){
				LL-=0.5*pow(modpars->Theta[ipar]*modpars->GSCALE,2);
			}
		}
	}
	/*
	if(LL>bestLL){
		bestLL=LL;
		CLog::Info("-------------------------------\n");
		CLog::Info("Theta=");
		for(ipar=0;ipar<NPars;ipar++){
			CLog::Info(to_string(modpars->Theta[ipar])+" ");
		}
		CLog::Info("\nY = ");
		for(iy=0;iy<NObs;iy++){
			CLog::Info(to_string(Y[iy])+"=?"+to_string(obsinfo->YExp[iy])+" ");
		}
		CLog::Info("\nbestLL="+to_string(bestLL)+"\n");
	}
	*/
	//Misc::Pause();
}

void CLLCalcSmooth::CalcLLPlusDerivatives(CModelParameters *modpars,double &LL,vector<double> &dLL_dTheta){
	unsigned int iy,ipar;
	double sigma2,delY;
	vector<vector<double>> dYdTheta;
	dYdTheta.resize(NObs);
	for(iy=0;iy<NObs;iy++){
		dYdTheta[iy].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			dYdTheta[iy][ipar]=0.0;
	}
	master->CalcAllYdYdTheta(modpars,Y,SigmaY_emulator,dYdTheta);
	LL=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		dLL_dTheta[ipar]=0.0;
	for(iy=0;iy<NObs;iy++){
		sigma2=SigmaY_emulator[iy]*SigmaY_emulator[iy]+obsinfo->SigmaExp[iy]*obsinfo->SigmaExp[iy];
		delY=Y[iy]-obsinfo->YExp[iy];
		LL-=0.5*delY*delY/sigma2;
		for(ipar=0;ipar<NPars;ipar++){
			dLL_dTheta[ipar]-=delY*dYdTheta[iy][ipar]/sigma2;
		}
	}
	for(ipar=0;ipar<NPars;ipar++){
		if(modpars->priorinfo->type[ipar]=="gaussian"){
			LL-=0.5*pow(modpars->Theta[ipar]*modpars->GSCALE,2);
			dLL_dTheta[ipar]-=modpars->Theta[ipar]*pow(modpars->GSCALE,2);
		}
	}
	/*
	if(LL>bestLL){
		bestLL=LL;
		CLog::Info("-------------------------------\n");
		CLog::Info("Theta=");
		for(ipar=0;ipar<NPars;ipar++){
			CLog::Info(to_string(modpars->Theta[ipar])+" ");
		}
		CLog::Info("\nY = ");
		for(iy=0;iy<NObs;iy++){
			CLog::Info(to_string(Y[iy])+"=?"+to_string(obsinfo->YExp[iy])+" ");
		}
		CLog::Info("\nbestLL="+to_string(bestLL)+"\n");
	}
	*/
	
	
}