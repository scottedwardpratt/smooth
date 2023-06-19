
#include<cstdlib>
#include <fstream>
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
    string s,s1;
    int i,Npars;
    double min,max,exp_val;
    vector<double> min_vec;
    vector<double> max_vec;
    //CPriorInfo* pInfo = new CPriorInfo("Info/mod_parametes_info.txt");
  //  CModelParameters modPar = CModelParameters(pInfo);
  //  Npars=modPar.NModelPars;

    FILE *fptr;
    fptr = fopen("obs.txt","w");

    FILE *fp;
    fp = fopen ("modelruns/run0/mod_parameters.txt","r");
    if (fp!=NULL)
    {
      for(int i = 1; i < 10; i ++)
      {
          fscanf(fp,"   %s %(s): x=%f, theta=%f",&s,&s,&min,&max);
          min_vec.push_back(min);
          max_vec.push_back(min);
      }
      for (int i = 0; i < 10; i++) {
        exp_val = min_vec[i]+max_vec[i];

        fprintf(fptr,"obs%d %s %10.3f %10.3f\n",i,s.c_str(),exp_val,0.001);
      }
      fclose (fp);
    }

    //cout <<
    fclose(fptr);
    return 0;


}
