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
    double min,max;
    FILE *fptr;
    fptr = fopen("obs.txt","w");

    FILE *fp;
    fp = fopen ("Info/mod_parameters_info.txt","r");
    if (fp!=NULL)
    {
      for(int i = 1; i < 10; i ++)
      {
          fscanf(fp,"par%d %s      %10.3f      %10.3f\n",&i,&s,&min,&max);
      }
      fclose (fp);
    }
    for(int i = 1; i < 10; i ++)
    {
        fprintf(fptr,"par%d %s %10.3f %10.3f\n",i,s.c_str(),min,max);
    }

    //cout <<
    fclose(fptr);
    return 0;


}
