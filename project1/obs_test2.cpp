
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
    string s;
    int i,Npars=9;
    double min,max,exp_val;
    vector<double> min_vec;
    vector<double> max_vec;
    //CPriorInfo* pInfo = new CPriorInfo("Info/mod_parametes_info.txt");
  //  CModelParameters modPar = CModelParameters(pInfo);
  //  Npars=modPar.NModelPars;

    FILE *fptr;
    //fptr = fopen("obs.txt","w");
    for(int ipars = 0; ipars < Npars; ipars++)
    {
      string dirname = "modelruns/run" + std::to_string(ipars);

    FILE *fp;
    fp = fopen ("/mod_parameters.txt","r");
    if (fp!=NULL)
    {
      for(int i = 1; i < 10; i ++)
      {
          fscanf(fp,"   %s : x=%lf, theta=%lf",&s,&s,&min,&max);
          min_vec.push_back(min);
          max_vec.push_back(min);
      }
      for (int i = 0; i < 10; i++) {
        FILE *fptr;

        string filename ="modelruns/run" + to_string(i)+ "/obs_test.txt";
        fptr = fopen(filename.c_str(),"a");
        exp_val = min_vec[i]+min_vec[i+1];

        fprintf(fptr,"obs%d %s %10.3f %10.3f\n",i,s.c_str(),exp_val,0.001);

      }
      fclose (fp);
    }
    //cout << min_vec[2] <<endl;
  }
    fclose(fptr);

    cout << min_vec[2] <<endl;
    return 0;

  }

  string s;
  int i;
  double min,max,exp_val;
  vector<double> min_vec;
  //vector<double> exp_val;
  FILE *fptr;
  fptr = fopen("obs.txt","w");

  FILE *fp;
  fp = fopen ("Info/mod_parameters_info.txt","r");
  if (fp!=NULL)
  {
    for(int i = 1; i < 10; i ++)
    {
        fscanf(fp,"par%d %s %lf %lf\n",&i,&s,&min,&max);
        min_vec.push_back(min);

        exp_val = min_vec[3];





    }
    fprintf(fptr,"obs%d %s %10.3f %10.3f\n",i,s.c_str(),exp_val,0.001);
    cout << min_vec[1] << endl;
    fclose (fp);
  }

  //cout <<
  fclose(fptr);
  return 0;
>>>>>>> cc96b5d (nsn)
