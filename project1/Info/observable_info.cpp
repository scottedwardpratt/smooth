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
    double exp_val, uncertainity;
    FILE *fptr;
    fptr = fopen("observable_info.txt","w");


    srand((unsigned)time(NULL));

    double step = 0.1;
    for(int i = 1; i < 10; i ++)
    {
        exp_val = (double)(rand()+ 0.5)/RAND_MAX;
        fprintf(fptr,"obs%d %10.3f %10.3f\n",i,exp_val,0.001);
    }
    fclose(fptr);


}
