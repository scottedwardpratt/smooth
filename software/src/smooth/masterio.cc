#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

void CSmoothMaster::ReadTrainingInfo(){
	if(SmoothEmulator_TrainingFormat == "SMOOTH"){
		traininginfo->ReadTrainingInfoSmoothFormat();
	}
	else if(SmoothEmulator_TrainingFormat == "SURMISE"){
		traininginfo->ReadTrainingInfoSurmiseFormat();
	}
}

void CSmoothMaster::ReadTestingInfo(){
   if(SmoothEmulator_TestingFormat == "SMOOTH"){
      testinginfo->ReadTestingInfoSmoothFormat();
   }
   else if(SmoothEmulator_TestingFormat == "SURMISE"){
      testinginfo->ReadTestingInfoSurmiseFormat();
   }
}

void CSmoothMaster::WriteSigmaLambda(string filename){
   FILE *fptr=fopen(filename.c_str(),"w");
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      fprintf(fptr,"%24s %10.3f %10.5f\n",
              observableinfo->observable_name[iY].c_str(),emulator[iY]->SigmaA,emulator[iY]->LAMBDA);
   }
   fclose(fptr);
}

void CSmoothMaster::WriteSigmaLambda(){
   string filename="smooth_data/output_stuff/sigmalambda.txt";
   WriteSigmaLambda(filename);
}

void CSmoothMaster::ReadSigmaLambda(){
   string filename="smooth_data/output_stuff/sigmalambda.txt";
   ReadSigmaLambda(filename);
}

void CSmoothMaster::ReadSigmaLambda(string filename){
   FILE *fptr=fopen(filename.c_str(),"r");
   char obsname[200];
   string obsstring;
   for(unsigned int iY=0;iY<observableinfo->NObservables;iY++){
      fscanf(fptr,"%s %lf %lf",obsname,&emulator[iY]->SigmaA,&emulator[iY]->LAMBDA);
      obsstring=obsname;
      if(obsstring!=observableinfo->observable_name[iY]){
         CLog::Fatal("Reading in obs name ("+obsstring+") from "+filename+" - does not match name in observable info ("+observableinfo->observable_name[iY]+")");
      }
   }
   fclose(fptr);
}
