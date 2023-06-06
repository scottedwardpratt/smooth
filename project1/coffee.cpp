#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include<vector>
using namespace std;
#include<cstdlib>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

class Run
{
public:
  static
  int count;
  int no;

  Run()
  {
    no = count;
    count++ ;
  }

  void MakeDir()
  {
    int stat;
    string path_s = "git/smooth/project1/modelruns/run%d", no;

    stat = _mkdir(path_s);

    if(stat == -1){
      cout << "Failed " << endl;
    }
  }

};

using namespace std;

int main()
{
  Run new_run;
  new_run = new Run();
  new_run.MakeDir();

  return 0;
}
