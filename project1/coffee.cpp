#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>

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

  return 0;
}
