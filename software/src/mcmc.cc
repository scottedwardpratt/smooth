#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

CMCMC::CMCMC(){
	//
}

CMCMC::CMCMC(CSmoothMaster *master_set){
	stepsize=0.02;
	master=master_set;
	parmap=master->parmap;
	randy=master->randy;
	NPars=master->NPars;
	
	ClearTrace();
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		trace[0].Theta[ipar]=0;
	}
	trace[0].TranslateTheta_to_X();
	llcalc=new CLLCalcSmooth(master);
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

void CMCMC::PerformMetropolisTrace(unsigned int Ntrace,unsigned int NSkip){
	unsigned int itrace,iskip,ipar;
	double oldLL,newLL,X;
	CModelParameters *oldmodpars,*newmodpars,*modpars;
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}

	for(itrace=0;itrace<Ntrace;itrace++){
		oldmodpars=&trace[trace.size()-1];
		newmodpars=new CModelParameters();
		llcalc->CalcLL(oldmodpars,oldLL);
		for(iskip=0;iskip<NSkip;iskip++){
			for(ipar=0;ipar<NPars;ipar++){
				newmodpars->Theta[ipar]=oldmodpars->Theta[ipar]+stepsize*randy->ran_gauss();
			}
			llcalc->CalcLL(newmodpars,newLL);
			if(newLL>oldLL){
				*oldmodpars=*newmodpars;
				oldLL=newLL;
			}
			if(newLL<oldLL){
				X=newLL-oldLL;
				if(X<-50)
					X=-50;
				if(exp(X)<randy->ran()){
					*oldmodpars=*newmodpars;
					oldLL=newLL;
				}
			}
		}
		modpars=new CModelParameters();
		*modpars=*oldmodpars;
		modpars->TranslateTheta_to_X();
		trace.push_back(*modpars);
		CLog::Info("finished for itrace"+to_string(itrace)+"\n");
		delete newmodpars;
	}
	
}