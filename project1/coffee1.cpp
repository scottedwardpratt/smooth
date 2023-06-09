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
  int iruns,Nruns=5;
  FILE *fptr;

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parametes_info.txt");
  CModelParameters modPar = CModelParameters(pInfo);
  modPar.TranslateX_to_Theta();
  modPar.TranslateTheta_to_x();


  for(int iruns = 0; iruns < Nruns; iruns++)
  {
    file >> iruns;
    string dirname = "modelruns/run" + std::to_string(iruns);
    string shellcommand = "mkdir -p "+dirname;
    system(shellcommand.c_str());

    string filename = "run" + to_string(iruns)+"mod_parametes.txt";
    fptr = fopen(filename.c_str(),"w");

    fprintf(fptr,"%11.4d\n",20);
    fclose(fptr);

  }

    modPar.Print();

  return 0;
}
