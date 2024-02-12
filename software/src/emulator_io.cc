#include "msu_smooth/emulator.h"
using namespace std;

using namespace NBandSmooth;
using namespace NMSUPratt;

void CSmoothEmulator::PrintA(vector<double> &Aprint){
	for(unsigned int ic=0;ic<smooth->NCoefficients;ic++){
		CLog::Info(to_string(ic)+": "+to_string(Aprint[ic])+"\n");
	}
}

void CSmoothEmulator::WriteCoefficients(){
	unsigned int isample,ic;
	FILE *fptr;
	string filename;
	string dirname=smoothmaster->CoefficientsDirName+"/"+observable_name;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	filename=dirname+"/meta.txt";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"# NPars  MaxRank NCoefficients\n");
	fprintf(fptr,"%u  %u  %u\n",NPars,smooth->MaxRank,smooth->NCoefficients);
	fclose(fptr);
	
	if(TuneChooseExact){
		filename=dirname+"AExact.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ic=0;ic<smooth->NCoefficients;ic++){
			fprintf(fptr,"%15.8e\n",AExact[ic]);
		}
		fclose(fptr);
	}
	else{
		for(isample=0;isample<NASample;isample++){
			filename=dirname+"/sample"+to_string(isample)+".txt";
			fptr=fopen(filename.c_str(),"w");
			for(ic=0;ic<smooth->NCoefficients;ic++){
				fprintf(fptr,"%15.8e\n",ASample[isample][ic]);
			}
			fclose(fptr);
		}
	}
}

void CSmoothEmulator::ReadCoefficients(){
	unsigned int isample,ic,NPars_test,MaxRank_test,NC_test;
	FILE *fptr;
	string filename;
	string dirname=smoothmaster->CoefficientsDirName+"/"+observable_name;
	string command="mkdir -p "+dirname;
	char dummy[100];
	system(command.c_str());
	filename=dirname+"/meta.txt";
	fptr=fopen(filename.c_str(),"r");
	fgets(dummy,100,fptr);
	fscanf(fptr,"%u  %u  %u\n",&NPars_test,&MaxRank_test,&NC_test);
	if(NPars_test!=NPars || MaxRank_test!=smooth->MaxRank || NC_test!=smooth->NCoefficients){
		CLog::Fatal("Mismatch in array sizes in ReadCoefficients");
	}
	fclose(fptr);
	
	if(TuneChooseExact){
		filename=dirname+"/AExact.txt";
		fptr=fopen(filename.c_str(),"r");
		for(ic=0;ic<smooth->NCoefficients;ic++){
			fscanf(fptr,"%lf\n",&ASample[isample][ic]);
		}
	}
	else{
		for(isample=0;isample<NASample;isample++){
			filename=dirname+"/sample"+to_string(isample)+".txt";
			fptr=fopen(filename.c_str(),"r");
			for(ic=0;ic<smooth->NCoefficients;ic++){
				fscanf(fptr,"%lf\n",&ASample[isample][ic]);
			}
			fclose(fptr);
		}
	}

}
