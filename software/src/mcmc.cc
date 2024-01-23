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
	trace_filename=parmap->getS("MCMC_TRACE_FILENAME","mcmc_traces/trace.txt");
	stepsize=parmap->getD("MCMC_METROPOLIS_STEPSIZE",0.05);
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

// erases trace, but keeps last point
void CMCMC::PruneTrace(){
	CModelParameters modpars;
	modpars=trace[trace.size()-1];
	trace.clear();
	trace.push_back(modpars);
	trace[0].TranslateTheta_to_X();
}

//void CMCMC::BurnInMetropolis(unsigned int Nburn){
//	for(unsigned int iburn=0;iburn<Nburn;iburn++){
//		Metropolis(NSkip);
//	}
//}

void CMCMC::PerformMetropolisTrace(unsigned int Ntrace,unsigned int NSkip){
	unsigned int nsuccess=0;
	unsigned int itrace,iskip,ipar,it0;
	double oldLL,newLL,X;
	CModelParameters *oldptr,*newptr;
	CModelParameters newmodpars;
	it0=trace.size();
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	trace.resize(trace.size()+Ntrace);
	for(itrace=it0;itrace<it0+Ntrace;itrace++){
		oldptr=&trace[trace.size()-1];
		newptr=&newmodpars;
		llcalc->CalcLL(oldptr,oldLL);
		for(iskip=0;iskip<NSkip;iskip++){
			for(ipar=0;ipar<NPars;ipar++){
				newptr->Theta[ipar]=oldptr->Theta[ipar]+stepsize*randy->ran_gauss();
			}
			llcalc->CalcLL(newptr,newLL);
			if(newLL>oldLL){
				oldptr->Copy(newptr);
				oldLL=newLL;
				nsuccess+=1;
			}
			if(newLL<oldLL){
				X=newLL-oldLL;
				if(X<-50)
					X=-50;
				if(exp(X)>randy->ran()){
					oldptr->Copy(newptr);
					oldLL=newLL;
					nsuccess+=1;
				}
			}
		}
		trace[itrace].TranslateTheta_to_X();
		trace[itrace].Copy(oldptr);
	}
	CLog::Info("Metropolis success percentage="+to_string(100.0*double(nsuccess)/(double(Ntrace)*double(NSkip)))+"\n");
}

void CMCMC::PerformLangevinTrace(unsigned int Ntrace,unsigned int NSkip){
	unsigned int itrace,iskip,ipar;
	bool inside;
	double LL,ss,sqstep,dLLdTheta2;
	vector<double> dLLdTheta,dTheta;
	CModelParameters *oldmodpars,*newmodpars,*modpars;
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	dLLdTheta.resize(NPars);
	dTheta.resize(NPars);
	
	for(itrace=0;itrace<Ntrace;itrace++){
		oldmodpars=&trace[trace.size()-1];
		newmodpars=new CModelParameters();
		for(iskip=0;iskip<NSkip;iskip++){
			//double oldLL=LL;
			llcalc->CalcLLPlusDerivatives(oldmodpars,LL,dLLdTheta);
			//double LLa,LLb,dLL;
			//llcalc->CalcLL(oldmodpars,LLa);
			/*
			if(LL<oldLL){
				printf("wrong way, oldLL=%g, LL=%g, dTheta=\n",oldLL,LL);
				for(ipar=0;ipar<NPars;ipar++){
					printf("%9.6f ",dTheta[ipar]);
				}
				printf("\n");
			}
			*/
			inside=true;
			dLLdTheta2=0.0;
			for(ipar=0;ipar<NPars;ipar++){
				dLLdTheta2+=dLLdTheta[ipar]*dLLdTheta[ipar];
			}
			ss=stepsize/sqrt(dLLdTheta2);
			sqstep=sqrt(2.0*ss);
			for(ipar=0;ipar<NPars;ipar++){
				dTheta[ipar]=ss*dLLdTheta[ipar]+sqstep*randy->ran_gauss();
				if(dTheta[ipar]!=dTheta[ipar]){
					oldmodpars->Print();
					CLog::Fatal("disaster in PerformLangevinTrace\n");
				}
				newmodpars->Theta[ipar]=oldmodpars->Theta[ipar]+dTheta[ipar];
				if(fabs(newmodpars->Theta[ipar])>1.0){
					inside=false;
					//ipar=NPars;
				}
			}
			
			if(inside){
				/*
				llcalc->CalcLL(newmodpars,LLb);
				newmodpars->Print();
				dLL=0.0;
				for(ipar=0;ipar<NPars;ipar++){
					dLL+=dLLdTheta[ipar]*dTheta[ipar];
					printf("%2u: dLLdTheta=%g, dTheta=%g\n",ipar,dLLdTheta[ipar],dTheta[ipar]);
				}
				printf("LLa=%g, LLb=%g, LL=%g\n",LLa,LLb,LL);
				printf("dLL=%g =? %g\n",dLL,LLb-LLa);
				Misc::Pause();
				*/
				*oldmodpars=*newmodpars;
			}
		}
		modpars=new CModelParameters();
		*modpars=*oldmodpars;
		modpars->TranslateTheta_to_X();
		trace.push_back(*modpars);
		//CLog::Info("finished for itrace"+to_string(itrace)+"\n");
		delete newmodpars;
	}
	
}

void CMCMC::WriteTrace(){
	unsigned int itrace,ipar,ntrace=trace.size();
	FILE *fptr;
	printf("writing, ntrace=%u\n",ntrace);
	fptr=fopen(trace_filename.c_str(),"w");
	for(itrace=0;itrace<ntrace;itrace++){
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%8.5f ",trace[itrace].Theta[ipar]);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
}

void CMCMC::ReadTrace(){
	unsigned int ipar;
	double thetaread;
	CModelParameters *modpars;
	FILE *fptr;
	fptr=fopen(trace_filename.c_str(),"r");
	while(!feof(fptr)){
		for(ipar=0;ipar<NPars;ipar++){
			fscanf(fptr,"%lf",&thetaread);
			if(ipar==0 && !feof(fptr)){
				modpars=new CModelParameters();
				modpars->TranslateTheta_to_X();
				trace.push_back(*modpars);
			}
			if(!feof(fptr)){
				modpars->Theta[ipar]=thetaread;
			}
		}
	}
	fclose(fptr);
	CLog::Info("read in "+to_string(trace.size())+" trace elements\n");
}
