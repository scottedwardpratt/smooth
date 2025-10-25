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

