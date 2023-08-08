#include "msu_smooth/modelparinfo.h"

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
			if(string(dummy2)!="uniform" && string(dummy2)!="gaussian"){
				CLog::Fatal("reading priorinfo: type="+string(dummy2)+". Must be uniform or gaussian.\n");
			}
			parname.push_back(string(dummy1));
			type.push_back(string(dummy2));
			xmin.push_back(minval);  // Gaussian type xmin refers to <x> and <xmax> is sigma
			xmax.push_back(maxval);
			name_map.insert(pair<string,int>(parname[NModelPars],NModelPars));
			NModelPars+=1;
		}
	}while(!feof(fptr));
	fclose(fptr);
}


int CPriorInfo::GetIPosition(string par_name){
	map<string,int>::iterator iter;
	pair<string,int> mpair;
	iter=name_map.find(par_name);
	if(iter==name_map.end()){
		CLog::Fatal("In CPriorInfo::GetIposition, cannot find parameter "+par_name+"\n");
	}
	return iter->second;
}

string CPriorInfo::GetName(int i){
	return parname[i];
}


CModelParameters::CModelParameters(CPriorInfo *priorinfo_set){
	priorinfo=priorinfo_set;
	NModelPars=priorinfo->NModelPars;
	X.resize(NModelPars);
	Theta.resize(NModelPars);
};

void CModelParameters::TranslateX_to_Theta(){
	//for min, max range
	double sigmax,xbar;
	int ipar;

	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			Theta[ipar]=-1+2*((X[ipar]-priorinfo->xmin[ipar])/(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]));
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			Theta[ipar]=(X[ipar]-xbar)/sigmax;
		}
		else{
			CLog::Fatal("Cannot translate X to Theta because type = "+priorinfo->type[ipar]+" is not recognized\n");
		}
	}
}

void CModelParameters::TranslateTheta_to_X(){
	double sigmax,xbar;
	int ipar;

	for(ipar=0;ipar<NModelPars;ipar++){
		if(priorinfo->type[ipar]=="uniform"){
			X[ipar]=priorinfo->xmin[ipar]+0.5*(1.0+Theta[ipar])*(priorinfo->xmax[ipar]-priorinfo->xmin[ipar]);
		}
		else if(priorinfo->type[ipar]=="gaussian"){
			xbar=priorinfo->xmin[ipar];
			sigmax=priorinfo->xmax[ipar];
			X[ipar]=xbar+sigmax*Theta[ipar];
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
		priorinfo->parname[ipar].c_str(),priorinfo->type[ipar].c_str(),X[ipar],Theta[ipar]);
		CLog::Info(message);
	}
}
