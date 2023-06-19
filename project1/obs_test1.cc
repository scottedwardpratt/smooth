#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include<vector>
using namespace std;
#include<cstdlib>
#include <fstream>





int main()
{
  string s;
  int i;
  double min,max,vals;
  vector<double> obs;
  //FILE *fptr;
  //fptr = fopen("obs.txt","w");

  FILE *fp;
  fp = fopen ("Info/mod_parameters_info.txt","r");
  if (fp!=NULL)
  {
    for(int i = 1; i < 10; i ++)
    {
      fscanf(fp,"par%d %s %lf %lf\n",&i,&s,&min,&max);
      while(min == 0.336) //loop forever
      {
          obs.push_back(min); //otherwise add it to the vector
      }


      //for(int i = 1; i<10; i= i+3 ){
      cout << obs[1]<< endl;
      //}
      //fprintf(fptr,"obs%d %s %10.3f %10.3f\n",i,s.c_str(),10.99,0.001);
    }
    fclose (fp);
  }

  //cout <<
  //fclose(fptr);
  return 0;


}
