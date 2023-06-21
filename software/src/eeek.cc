#include "msu_smooth/real.h"
#include "msu_smooth/smooth.h"
#include "msu_commonutils/randy.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"

using namespace std;

//Creating a real fucntion for Project1.


CReal::CReal(){
//	cout << "CReal object created" << endl;
}

void CReal::CalcY(vector<double> &theta,double &Y,double &SigmaY)
{
	cout << "dummy function -- should not be hear" << endl;
	Y=SigmaY=0.0;
}

CReal_EEEK::CReal_EEEK(unsigned int NPars_Set)
{
	NPars=NPars_Set;
	LAMBDA=2;
}


void CReal_EEEK::CalcY(vector<double> &theta,double &Y,double &SigmaY)
{
	Y=smooth->CalcY(A,LAMBDA,theta);
	SigmaY=1.0;
}

void CReal::CalcYTrain(vector<double> &YTrain,vector<double> &SigmaYTrain, int NTrainingPts, vector<vector<double>> ThetaTrain){
//	cout << "NTrainingPts" << NTrainingPts << endl;
	//NtrainingPts is 4
	unsigned int itrain;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		CalcY(ThetaTrain[itrain],YTrain[itrain],SigmaYTrain[itrain]);
	}
}
