#include "msu_smooth/master.h"
#include "msu_smooth/mcmc.h"
#include "msu_smoothutils/misc.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CMCMC::CMCMC(){
	//
}

CMCMC::CMCMC(CSmoothMaster *master_set){
	master=master_set;
	parmap=master->parmap;
	randy=master->randy;
	NPars=master->NPars;
	trace_filename=parmap->getS("MCMC_TRACE_FILENAME","mcmc_trace/trace.txt");
	OPTIMIZESTEPS=parmap->getB("MCMC_OPTIMIZESTEPS",false);
	langevin=parmap->getB("MCMC_LANGEVIN",false);
	if(langevin)
		stepsize=parmap->getD("MCMC_LANGEVIN_STEPSIZE",0.01);
	else
		stepsize=parmap->getD("MCMC_METROPOLIS_STEPSIZE",0.05);
	ClearTrace();
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		trace[0].Theta[ipar]=0.0;
	}
	trace[0].TranslateTheta_to_X();
	llcalc=new CLLCalcSmooth(master);
	
	OPTIMIZESTEPS=false;
	dThetadTheta.resize(NPars,NPars);
	dTdTEigenVecs.resize(NPars,NPars);
	dTdTEigenVals.resize(NPars);
	stepvec.resize(NPars);
	stepvecprime.resize(NPars);
	dTdTEigenVecs.setZero();
	for(unsigned int ipar=0;ipar<NPars;ipar++){
		dTdTEigenVecs(ipar,ipar)=1.0;
		dTdTEigenVals(ipar)=stepsize;
	}
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

void CMCMC::PerformTrace(unsigned int Ntrace,unsigned int Nskip){
	if(langevin)
		PerformLangevinTrace(Ntrace,Nskip);
	else
		PerformMetropolisTrace(Ntrace,Nskip);
}

void CMCMC::PerformMetropolisTrace(unsigned int Ntrace,unsigned int Nskip){
	unsigned long long int nsuccess=0;
	unsigned int itrace,iskip,ipar,it0;
	double oldLL,newLL,X,bestLL;
	CModelParameters *oldptr,*newptr;
	CModelParameters oldmodpars;
	it0=trace.size();
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	OPTIMIZESTEPS=false;
	
	it0=trace.size();
	trace.resize(trace.size()+Ntrace);
	if(it0==0){
		oldmodpars.X.resize(NPars);
		oldmodpars.priorinfo=master->priorinfo;
		oldptr=&oldmodpars;
		for(ipar=0;ipar<NPars;ipar++){
			oldmodpars.Theta[ipar]=0.0;
		}
		oldptr->TranslateTheta_to_X();
	}
	else{
		oldmodpars.Copy(&trace[it0-1]);
	}
	oldptr=&oldmodpars;
	oldptr->TranslateTheta_to_X();
	
	llcalc->CalcLL(oldptr,oldLL);
	//oldLL=0.0;
	//for(ipar=0;ipar<NPars;ipar++)
	//	oldLL-=pow(oldptr->Theta[ipar]-0.2,2);
	bestLL=oldLL;
	CLog::Info("At beginning of Trace, LL="+to_string(oldLL)+"\n");
	
	for(itrace=it0;itrace<it0+Ntrace;itrace++){
		for(iskip=0;iskip<Nskip;iskip++){
			newptr=&trace[itrace];
				
			if(OPTIMIZESTEPS){
				for(ipar=0;ipar<NPars;ipar++){
					stepvecprime[ipar]=stepsize*randy->ran_gauss();
				}
				stepvec=dTdTEigenVecs*stepvecprime;
				for(ipar=0;ipar<NPars;ipar++){
					newptr->Theta[ipar]=oldptr->Theta[ipar]+real(stepvec(ipar));
				}
			}
			else{
				for(ipar=0;ipar<NPars;ipar++){
					newptr->Theta[ipar]=oldptr->Theta[ipar]+stepsize*randy->ran_gauss();
				}
			}
			
			llcalc->CalcLL(newptr,newLL);
			//newLL=0.0;
			//for(ipar=0;ipar<NPars;ipar++)
			//	newLL-=10*pow(newptr->Theta[ipar]-0.2,2);
			
			if(newLL>bestLL){
				bestLL=newLL;
			}
			
			if(newLL>oldLL){
				oldptr->Copy(newptr);
				oldLL=newLL;
				nsuccess+=1;
			}
			else{
				X=newLL-oldLL;
				if(X>-30.0){
					if(exp(X)>randy->ran()){
						oldptr->Copy(newptr);
						oldLL=newLL;
						nsuccess+=1;
					}
				}
			}
		}
		oldptr->TranslateTheta_to_X();
		trace[itrace].Copy(oldptr);
		if(Ntrace>10 && ((itrace+1)*10)%Ntrace==0){
			CLog::Info("finished "+to_string(lrint(100*double(itrace+1)/double(Ntrace)))+"%\n");
		}
	}
	CLog::Info("At end of trace, best LL="+to_string(bestLL)+"\n");
	CLog::Info("Metropolis success percentage="+to_string(100.0*double(nsuccess)/(double(Ntrace)*double(Nskip)))+"\n");
}

void CMCMC::PerformLangevinTrace(unsigned int Ntrace,unsigned int Nskip){
	unsigned int itrace,iskip,ipar;
	bool inside;
	double LL,ss,sqstep,dLLdTheta2,bestLL=-1.0E99;
	vector<double> dLLdTheta,dTheta;
	CModelParameters *oldptr,*newmodpars,*modpars;
	if(trace.size()==0){
		CLog::Fatal("Inside MCMC::PerformMetropolis, no initial point!\n");
	}
	dLLdTheta.resize(NPars);
	dTheta.resize(NPars);
	llcalc->CalcLL(&trace[trace.size()-1],bestLL);
	CLog::Info("At beginning of trace LL="+to_string(bestLL)+"\n");
	
	for(itrace=0;itrace<Ntrace;itrace++){
		oldptr=&trace[trace.size()-1];
		newmodpars=new CModelParameters();
		for(iskip=0;iskip<Nskip;iskip++){
			llcalc->CalcLLPlusDerivatives(oldptr,LL,dLLdTheta);
			if(LL>bestLL){
				bestLL=LL;
			}
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
					oldptr->Print();
					CLog::Fatal("disaster in PerformLangevinTrace\n");
				}
				newmodpars->Theta[ipar]=oldptr->Theta[ipar]+dTheta[ipar];
				if(fabs(newmodpars->Theta[ipar])>1.0){
					inside=false;
					//ipar=NPars;
				}
			}
			
			if(inside){
				*oldptr=*newmodpars;
			}
		}
		modpars=new CModelParameters();
		*modpars=*oldptr;
		modpars->TranslateTheta_to_X();
		trace.push_back(*modpars);
		//CLog::Info("finished for itrace"+to_string(itrace)+"\n");
		delete newmodpars;
		if(Ntrace>10 && ((itrace+1)*10)%Ntrace==0){
			CLog::Info("finished "+to_string(lrint(100*double(itrace+1)/double(Ntrace)))+"%\n");
		}
	}
	CLog::Info("At end of trace: best LL="+to_string(bestLL)+"\n");
	
}

void CMCMC::WriteTrace(){
	unsigned int itrace,ipar,ntrace=trace.size();
	FILE *fptr;
	CLog::Info("writing, ntrace = "+to_string(ntrace)+"\n");
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
