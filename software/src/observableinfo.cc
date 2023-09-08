#include "msu_smooth/priorinfo.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <cstdio>
#include "msu_commonutils/log.h"
#include "msu_smooth/observableinfo.h"

using namespace std;

CObservableInfo::CObservableInfo(string filename){
	ReadObservableInfo(filename);
}

int CObservableInfo::GetIPosition(string obsname){
	map<string,int>::iterator iter;
	pair<string,int> mpair;
	iter=name_map.find(obsname);
	if(iter==name_map.end()){
		CLog::Fatal("In CObservableInfo::GetIposition, cannot find observable "+obsname+"\n");
	}
	return iter->second;
} 

string CObservableInfo::GetName(int i){
	return observable_name[i];
}

void CObservableInfo::ReadObservableInfo(string filename){
	char dummy1[120];
	double sig0;
	name_map.clear();
	observable_name.clear();
	NObservables=0;
	FILE *fptr=fopen(filename.c_str(),"r");
	do{
		fscanf(fptr,"%s %lf",dummy1,&sig0);
		if(!feof(fptr)){
			observable_name.push_back(string(dummy1));
			SigmaA0.push_back(sig0);
			name_map.insert(pair<string,int>(observable_name[NObservables],NObservables));
			NObservables+=1;
		}
	}while(!feof(fptr));
	
	printf("Nobs=%d\n",NObservables);
	for(int iy=0;iy<NObservables;iy++){
		printf("%s  %g\n",observable_name[iy].c_str(),SigmaA0[iy]);
	}
	//exit(1);
	fclose(fptr);
}

void CObservableInfo::ReadExperimentalInfo(string filename){
	char dummy1[120];
	double sig0,Y0;
	int NObsRead=0;
	int iY;
	YExp.resize(NObservables);
	SigmaExp.resize(NObservables);
	FILE *fptr=fopen(filename.c_str(),"r");
	do{
		fscanf(fptr,"%s %lf %lf",dummy1,&Y0,&sig0);
		if(!feof(fptr)){
			iY=GetIPosition(string(dummy1));
			YExp[iY]=Y0;
			SigmaExp[iY]=sig0;
			NObsRead+=1;
		}
	}while(!feof(fptr));
	fclose(fptr);
	if(NObsRead!=NObservables){
		CLog::Fatal("NObservables="+to_string(NObservables)+", N of Experimental Observables="+to_string(NObsRead)+"\n");
	}

}

void CObservableInfo::PrintInfo(){
	CLog::Info("Observable    \n");
	for(int i=0;i<NObservables;i++){
		CLog::Info(observable_name[i]+"\n");
	}
}