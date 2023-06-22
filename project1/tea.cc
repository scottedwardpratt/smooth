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
  std::cout << "Number of lines in text file: " << Nobs<< endl;
  return Nobs - 1;

}

int main()
{

  ifstream file;
  int Npars;
  int No_of_obs = noObs();
  FILE *fptr;
  unsigned int NPars;
  double y,yreal,accuracy,average_accuracy=0.0,average_expected_accuracy=0.0,sigmay2,ybar,y2bar,SigmaYreal;
	unsigned int isample,itest,ntest=10,ipar,ireal,nreal=1;
	vector<double> Theta;

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parameters_info.txt");
  CModelParameters modPar = CModelParameters(pInfo);



  CparameterMap *parmap = new CparameterMap();
  //double YExp,SigmaYExp,SigmaYReal;
  vector<vector<double>> ThetaTest;

  Npars=modPar.NModelPars;
  cout << "NUMBER OF parameter IS " << Npars <<endl;

  // This plays the role of the "real" model
  CReal_EEEK *real;

  //parmap->ReadParsFromFile(parfilename);
  parmap->set("SmoothEmulator_NPars",Npars);

  CSmoothEmulator emulator(parmap);
  emulator.randy->reset(-time(NULL));
  NPars = emulator.NPars = Npars;
  emulator.InitTrainingPtsArrays(Npars);

  //cout << "NUMBER OF parameter IS " << NPars <<endl;

  Theta.resize(NPars);

  // Real function
  real=new CReal_EEEK(emulator.NPars,emulator.smooth->MaxRank,emulator.randy);
  real->LAMBDA=emulator.LAMBDA;
  emulator.real=real;
  real->RandomizeA(100.0);


  for (size_t i_runs = 0; i_runs < Npars; i_runs++) {
    float obs;
    vector<float> obs_vals;

    FILE *fptr;
    string obs_dir = "modelruns/run" + to_string(i_runs) + "/obs.txt";

    fptr = fopen(obs_dir.c_str(), "r");

    for(int i = 0; i < No_of_obs; i++)
    {

      string s;
      fscanf(fptr,"%s %f\n",&s,&obs);
      cout << "hello" << endl;
      emulator.ThetaTrain[i_runs][i] = obs;
    }
    fclose (fptr);
  }

  real->RandomizeA(100.0);
  //real->A[0]=0.0;
  emulator.CalcYTrainFromThetaTrain();
  emulator.GenerateASamples();
  //  This is just for testing, to make sure that the emulator exactly reproduces training points
  for(int itest=0;itest<emulator.NTrainingPts;itest++){
    yreal=emulator.YTrain[itest];
    emulator.real->CalcY(emulator.ASample[1], emulator.ThetaTrain[itest], y, SigmaYreal);
    if(fabs(yreal-y)>0.001){
      printf("Emulator fails!\n");
      printf("%u, %8.4f: %g =? %g\n",itest,emulator.ThetaTrain[itest][0],y,yreal);
    }
  }
  //cout << pInfo <<endl;


  return 0;
}
