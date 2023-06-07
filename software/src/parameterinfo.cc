#include "msu_smooth/parameterinfo.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/log.h"

using namespace std;

CPriorInfo::CPriorInfo(string parinfo_filename_set){
	parinfo_filename=parinfo_filename_set;
	double minval,maxval;
	char dummy1[120],dummy2[40];
	
	FILE *fptr;
	fptr=fopen(parinfo_filename.c_str(), "r");
	NModelPars=0;
	do{
		fscanf(fptr, "%s %s %lf %lf",dummy1,dummy2,&minval,&maxval);
		if(!feof(fptr)){
			if(string(dummy2)!="linear" && string(dummy2)!="gaussian"){
				CLog::Fatal("reading priorinfo: type="+string(dummy2)+". Must be linear or gaussian.\n");
			}
			parname.push_back(string(dummy1));
			type.push_back(string(dummy2));
			xmin.push_back(minval);  // Gaussian type xmin refers to <x> and <xmax> is sigma
			xmax.push_back(maxval);
			NModelPars+=1;
		}
	}while(!feof(fptr));
	fclose(fptr);
}

CModelParameters::CModelParameters(CPriorInfo *priorinfo_set){
	priorinfo=priorinfo_set;
	NModelPars=priorinfo->NModelPars;
	x.resize(NModelPars);
	theta.resize(NModelPars);
};

void CModelParameters::TranslateX_to_Theta(){
	//for min, max range
	double sigmax,xbar;
	int ipar;
	
	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="linear"){
			theta[ipar]=-1+2*((x[ipar]-priorinfo->xmin[ipar])/(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]));
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			theta[ipar]=(x[ipar]-xbar)/sigmax;
		}
		else{
			CLog::Fatal("Cannot translate X to Theta because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}

}

void CModelParameters::TranslateTheta_to_x(){
	double sigmax,xbar;
	int ipar;
	
	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="linear"){
			x[ipar]=priorinfo->xmin[ipar]+0.5*(1.0+theta[ipar])*(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]);
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			x[ipar]=xbar+sigmax*theta[ipar];
		}
		else{
			CLog::Fatal("Cannot translate Theta to X because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}
	
}

void CModelParameters::Print(){
	char message[200];
	int ipar;
	for(ipar=0;ipar<NModelPars;ipar++){
		snprintf(message,200,"   %.24s (%.8s): x=%11.4e, theta=%11.4e\n",
		priorinfo->parname[ipar].c_str(),priorinfo->type[ipar].c_str(),x[ipar],theta[ipar]);
		CLog::Info(message);
	}
}
