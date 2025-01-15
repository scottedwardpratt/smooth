#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

void CSmoothMaster::ReadTrainingInfo(){
	if(SmoothEmulator_TrainingFormat == "training_format_smooth"){
		traininginfo->ReadTrainingInfoSmoothFormat();
	}
	else if(SmoothEmulator_TrainingFormat == "training_format_surmise"){
		traininginfo->ReadTrainingInfoSurmiseFormat();
	}
}

