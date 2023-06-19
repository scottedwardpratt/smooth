#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>

using namespace std;

int main()
{
  char s;
  float val,theta;
  vector<float> par_vals;
  vector<float> obs;

  for (size_t i = 0; i < 10; i++) {

  FILE *fptr;
  string obs_dir = "modelruns/run" + to_string(i) + "/obs.txt";

  FILE *fp;
  string dirname = "modelruns/run" + to_string(i) + "/mod_parameters.txt";
  fp = fopen(dirname.c_str(),"r");

    for(int i = 1; i < 10; i++)
    {
      fscanf(fp,"%s %f\n",&s,&val);
      par_vals.push_back(val);
    }


    for (size_t i = 0; i < 9; i++)
    {
      obs.push_back(par_vals[i] + par_vals[i+1]);
    }

    for (size_t i = 0; i < 3; i++)
    {
      obs.push_back(par_vals[i] * par_vals[i + 3]);
    }


    int obs_length = obs.size();

    for (size_t i = 0; i < obs_length; i++) {
      fptr = fopen(obs_dir.c_str(),"a");
      fprintf(fptr,"obs%ld %f",i, obs[i]);
      fclose (fptr);
    }
    fclose (fp);
  }
  return 0;


}
