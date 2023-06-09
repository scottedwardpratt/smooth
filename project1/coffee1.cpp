#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include "msu_smooth/parameterinfo.h"

using namespace std;

int main()
{
  ifstream file;
  int Npars;
  FILE *fptr;

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parametes_info.txt");
  CModelParameters modPar = CModelParameters(pInfo);
  modPar.TranslateX_to_Theta();
  modPar.TranslateTheta_to_x();

  Npars=modPar.NModelPars;


  for(int ipars = 0; ipars < Npars; ipars++)
  {
    file >> ipars;
    string dirname = "modelruns/run" + std::to_string(ipars);
    string shellcommand = "mkdir -p "+dirname;
    system(shellcommand.c_str());

    for(int i =0; i < Npars; i++)
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
