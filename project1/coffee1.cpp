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

int main()
{
  ifstream file;
  int Npars;
  FILE *fptr;
  CparameterMap *parmap=new CparameterMap();


  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parametes_info.txt");



  CModelParameters modPar = CModelParameters(pInfo);
  modPar.TranslateX_to_Theta();
  modPar.TranslateTheta_to_x();

  Npars=modPar.NModelPars;


  for(itest=0;itest<ntest;itest++){
    for(ipar=0;ipar<Npars;ipar++){
      ThetaTest[itest][ipar]=1.0-2.0*emulator.randy->ran();
    }
  }
  ThetaTest.resize(ntest);
  for(itest=0;itest<ntest;itest++){
    ThetaTest[itest].resize(NPars);
  }



  for(int ipars = 0; ipars < Npars+1; ipars++)
  {
    file >> ipars;
    string dirname = "modelruns/run" + std::to_string(ipars);
    string shellcommand = "mkdir -p "+dirname;
    system(shellcommand.c_str());

    for(int i =0; i < Npars+1; i++)
    {
      string filename ="modelruns/run" + to_string(ipars)+ "/mod_parametes.txt";
      fptr = fopen(filename.c_str(),"w");
      fprintf(fptr,"%11.4d\n",20);
      fclose(fptr);
    }

  }

    cout << modPar.NModelPars <<endl;
    modPar.Print();

  return 0;
}
