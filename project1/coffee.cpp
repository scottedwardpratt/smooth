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
  int runNo;
  FILE *fptr;
  ifstream file;

  file.open("log_file.txt");
  file >> runNo;

  string dirname = "modelruns/run" + std::to_string(runNo);
  string shellcommand = "mkdir -p "+dirname;
  system(shellcommand.c_str());

  fptr = fopen("log_file.txt","w");
  fprintf(fptr,"%d",runNo + 1);

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parametes_info.txt");

  CModelParameters modPar = CModelParameters(pInfo);
  modPar.TranslateX_to_Theta();
  modPar.TranslateTheta_to_x();

  modPar.Print();



  return 0;
}
