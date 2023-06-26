#include "msu_smooth/priorinfo.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/log.h"
#include "msu_smooth/observableinfo.h"

using namespace std;

CObservableInfo::CObservableInfo(string obs_inf_filename_set){
	observable_info_filename=obs_inf_filename_set;
	ReadObservableInfo(observable_info_filename);
}


int CObservableInfo::GetIPosition(string obsname){
	map<string,int>::iterator iter;
	pair<string,int> mpair;
	iter=name_map.find(obsname);
	if(iter==name_map.end()){
		CLog::Fatal("In CObservableInfo::GetIposition, cannot find observable"+obsname+"\n");
	}
	return iter->second;
} 

string CObservableInfo::GetName(int i){
	return observable_name[i];
}

void CObservableInfo::ReadObservableInfo(string observable_info_filename){
	char dummy1[120],dummy2[120];
	double sig0;
	name_map.clear();
	observable_name.clear();
	unit.clear();
	NObservables=0;
	FILE *fptr=fopen(observable_info_filename.c_str(),"r");
	do{
		fscanf(fptr,"%s %s %lf",dummy1,dummy2,&sig0);
		if(!feof(fptr)){
			observable_name.push_back(string(dummy1));
			unit.push_back(string(dummy2));
			SigmaA0.push_back(sig0);
			name_map.insert(pair<string,int>(observable_name[NObservables],NObservables));
			NObservables+=1;
		}
	}while(!feof(fptr));
	fclose(fptr);
}

void CObservableInfo::PrintInfo(){
	CLog::Info("Observable         Unit\n");
	for(int i=0;i<NObservables;i++){
		CLog::Info(observable_name[i]+"    "+unit[i]+"\n");
	}
}