#include <iostream>
using namespace std;
#include<cstdlib>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


int main()
{
    string dirname = "modelruns/run0";
    string shellcommand = "mkdir -p "+dirname;
    system(shellcommand.c_str());
    return 0;
}
