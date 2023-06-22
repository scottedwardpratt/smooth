#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include "msu_smooth/parameterinfo.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/scorecard.h"



using namespace std;

int noObs(){

  int Nobs = 0;
  std::string line;
  std::ifstream myfile("modelruns/run0/obs.txt");
  while (std::getline(myfile, line))
  ++Nobs;
  std::cout << "Number of lines in text file: " << Nobs - 1<< endl;
  return Nobs - 1;

}

int main()
{

  ifstream file;
  int Npars;
  int No_of_obs = noObs();
  int NPars;
  double y,yreal,accuracy,average_accuracy=0.0,average_expected_accuracy=0.0,sigmay2,ybar,y2bar,SigmaYreal;
	unsigned int isample,itest,ntest=25,ipar,ireal,nreal=1;
	vector<double> Theta;

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parameters_info.txt");
  CModelParameters modPar = CModelParameters(pInfo);

  CparameterMap *parmap = new CparameterMap();

  string parfilename = "parameters.txt";
  parmap->ReadParsFromFile(parfilename);
  //double YExp,SigmaYExp,SigmaYReal;
  vector<vector<double>> ThetaTest;

  NPars = Npars = modPar.NModelPars;
  cout << "NUMBER OF parameter IS " << Npars <<endl;

  // This plays the role of the "real" model
  CReal_EEEK *real;

  //parmap->ReadParsFromFile(parfilename);
  parmap->set("SmoothEmulator_NPars",Npars);

  CSmoothEmulator emulator(parmap);
  emulator.randy->reset(-time(NULL));
  emulator.NPars = No_of_obs;
  emulator.LAMBDA = 1;
  emulator.InitTrainingPtsArrays(NPars);
  Theta.resize(emulator.NPars);

  //cout << "NUMBER OF parameter IS " << NPars <<endl;

  real=new CReal_EEEK(emulator.NPars,emulator.smooth->MaxRank,emulator.randy);
  real->LAMBDA=emulator.LAMBDA;
  emulator.real=real;

  float obs;
  string s;
  FILE *fptr;
  string obs_dir;

  for (size_t i_runs = 0; i_runs < Npars; i_runs++) {
    obs_dir = "modelruns/run" + to_string(i_runs) + "/obs.txt";
    fptr = fopen(obs_dir.c_str(), "r");
    for(int i = 0; i < No_of_obs; i++)
    {
      fscanf(fptr,"%s %f\n",&s,&obs);
      emulator.ThetaTrain[i_runs][i] = obs;
    }
    fclose (fptr);
  }

  for(int ireal=0;ireal<nreal;ireal++){
    printf("------ ireal=%d -----\n",ireal);
    accuracy=0.0;
    real->RandomizeA(100.0);
    //real->A[0]=0.0;
    //  This is just for testing, to make sure that the emulator exactly reproduces training points
    emulator.CalcYTrainFromThetaTrain();
    emulator.GenerateASamples();

    for(int itest=0;itest < emulator.NTrainingPts;itest++){
      yreal=emulator.YTrain[itest];
      y=emulator.smooth->CalcY(emulator.ASample[1],emulator.LAMBDA,emulator.ThetaTrain[itest]);
      if(fabs(yreal-y)>0.001){
        printf("Emulator fails!\n");
        printf("%u, %8.4f: %g =? %g\n",itest,emulator.ThetaTrain[itest][0],y,yreal);
      }
    }

    accuracy=sigmay2=0.0;
		for(itest=0;itest<ntest;itest++){
			for(ipar=0;ipar<emulator.NPars;ipar++){
				if(emulator.NPars==1){
					Theta[ipar]=-1.0+(2.0/double(ntest))*(0.5+itest);
				}
				else{
					Theta[ipar]=1.0-2.0*emulator.randy->ran();
				}
			}
			real->CalcY(Theta,yreal,SigmaYreal);
			ybar=y2bar=0.0;
			for(isample=0;isample<emulator.NASample;isample++){
				y=emulator.smooth->CalcY(emulator.ASample[isample],emulator.LAMBDA,Theta);
				y2bar+=y*y;
				ybar+=y;
			}
			y2bar=y2bar/double(emulator.NASample);
			ybar=ybar/double(emulator.NASample);
			printf("ybar=%8.4f =? %8.4f\n",ybar,yreal);

			accuracy+=(ybar-yreal)*(ybar-yreal);
			sigmay2+=y2bar-ybar*ybar;
		}
		accuracy=accuracy/ntest;
		sigmay2=sigmay2/double(ntest-1);
		accuracy=sqrt(accuracy);
		sigmay2=sqrt(sigmay2);
		printf("accuracy=%7.3f, expected accuracy=%7.3f\n",accuracy,sigmay2);
		average_accuracy+=accuracy;
		average_expected_accuracy+=sigmay2;

	}
	average_accuracy=average_accuracy/double(nreal);
	average_expected_accuracy=average_expected_accuracy/double(nreal);

	printf("average <accuracy>=%g, average <expected accuracy>=%g\n",average_accuracy,average_expected_accuracy);

	return 0;
  }
