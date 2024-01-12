#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

CMCMC::CMCMC(){
	//
}

CMCMC::CMCMC(CSmoothMaster *master_set){
	master=master_set;
	parmap=master->parmap;
	randy=master->randy;
	NPars=master->NPars;
	ClearTrace();
	CModelParameters *modpars=new CModelParameters();
	trace.push_back(*modpars);
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		trace[0].Theta[ipar]=0;
	}
	trace[0].TranslateTheta_to_X();
}

void CMCMC::ClearTrace(){
	CModelParameters *modpars=new CModelParameters();
	trace.clear();
	trace.push_back(*modpars);
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		trace[0].Theta[ipar]=0;
	}
	trace[0].TranslateTheta_to_X();
}

void CMCMC::ClearBurnTrace(){
	burntrace.clear();
	CModelParameters *modpars=new CModelParameters();
	burntrace.push_back(*modpars);
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		burntrace[0].Theta[ipar]=0;
	}
	burntrace[0].TranslateTheta_to_X();
}

void CMCMC::ClearTrace(CModelParameters *modpars){
	trace.clear();
	trace.push_back(*modpars);
	trace[0].TranslateTheta_to_X();
}

void CMCMC::ClearBurnTrace(CModelParameters *modpars){
	burntrace.clear();
	burntrace.push_back(*modpars);
	burntrace[0].TranslateTheta_to_X();
}

//void CMCMC::BurnInMetropolis(unsigned int Nburn){
//	for(unsigned int iburn=0;iburn<Nburn;iburn++){
//		Metropolis(NSkip);
//	}
//}

void CMCMC::PerformMetropolisTrace(unsigned int Ntrace){
	unsigned int itrace,iskip,ipar;
	double oldLL,newLL;
	CModelParameters *oldmodpars,*newmodpars,*modpars;
	if(trace.size()<=0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	oldmodpars=&trace[trace.size()-1];
	llcalc->CalcLL(oldmodpars,oldLL);
	newmodpars=new CModelParameters();
	for(itrace=0;itrace<Ntrace;itrace++){
		for(iskip=0;iskip<NSkip;iskip++){
			for(ipar=0;ipar<NPars;ipar++){
				newmodpars->Theta[ipar]=oldmodpars->Theta[ipar]+stepsize*randy->ran_gauss();
			}
			llcalc->CalcLL(newmodpars,newLL);
			if(newLL>oldLL || exp(newLL-oldLL)>randy->ran()){
				*oldmodpars=*newmodpars;
			}
		}
		modpars=new CModelParameters();
		*modpars=*oldmodpars;
		modpars->TranslateTheta_to_X();
		trace.push_back(*modpars);
	}	
}