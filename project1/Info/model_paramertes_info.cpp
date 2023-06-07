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
    string parname, min_val, max_val;
    double min,max;
    FILE *fptr;
    fptr = fopen("mod_parametes_info.txt","w");


    srand((unsigned)time(NULL));

    double step = 0.1;
    for(int i = 1; i < 10; i ++)
    {
        min = i * step;
        max = min + i;
        fprintf(fptr,"par%d %10.3f %10.3f\n",i,min,max);
    }
    fclose(fptr);


}
