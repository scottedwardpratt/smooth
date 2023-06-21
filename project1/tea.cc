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
  std::cout << "Number of lines in text file: " << Nobs;
  return Nobs;

}

int main()
{

  ifstream file;
  int Npars;
  int No_of_obs = noObs();
  FILE *fptr;
  unsigned int NPars;

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parameters_info.txt");
  CModelParameters modPar = CModelParameters(pInfo);



  CparameterMap *parmap = new CparameterMap();
  //double YExp,SigmaYExp,SigmaYReal;
  vector<vector<double>> ThetaTest;
  vector<double> Theta;

  Npars=modPar.NModelPars;
  cout << "NUMBER OF parameter IS " << Npars <<endl;




  for (size_t i = 0; i < Npars; i++) {
    float obs;
    vector<float> obs_vals;

    FILE *fptr;
    string obs_dir = "modelruns/run" + to_string(i) + "/obs.txt";



    for(int i = 1; i < No_of_obs; i++)
    {
      string s;

      fscanf(fptr,"%s %f\n",&s,&obs);
      obs_vals.push_back(obs);
    }

    fclose (fptr);


  }



  // This plays the role of the "real" model
  CReal_EEEK *real;

  //parmap->ReadParsFromFile(parfilename);
  parmap->set("SmoothEmulator_NPars",Npars);

  CSmoothEmulator emulator(parmap);
  emulator.randy->reset(-time(NULL));
  NPars = emulator.NPars = Npars;

  //cout << "NUMBER OF parameter IS " << NPars <<endl;

  Theta.resize(NPars);

  // Real function
  real=new CReal_EEEK(emulator.NPars,emulator.smooth->MaxRank,emulator.randy);
  real->LAMBDA=emulator.LAMBDA;
  emulator.real=real;
  real->RandomizeA(100.0);








  //cout << pInfo <<endl;


  return 0;
}
