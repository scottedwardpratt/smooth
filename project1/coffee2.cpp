#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include "msu_smooth/parameterinfo.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/scorecard.h"





using namespace std;

int main()
{
  string removecommand = "rm -r modelruns/run*";
  system(removecommand.c_str());

  ifstream file;
  int Npars;
  FILE *fptr;

  CPriorInfo* pInfo = new CPriorInfo("Info/mod_parameters_info.txt");
  CModelParameters modPar = CModelParameters(pInfo);



  CparameterMap *parmap = new CparameterMap();
	//double YExp,SigmaYExp,SigmaYReal;
	unsigned int itest,ntest=1000,ipar,NPars;
	vector<vector<double>> ThetaTest;
	vector<double> Theta;

  Npars=modPar.NModelPars;
  cout << "NUMBER OF parameter IS " << Npars <<endl;


	// This plays the role of the "real" model
	CReal_Taylor *real;

	//parmap->ReadParsFromFile(parfilename);
  parmap->set("SmoothEmulator_NPars",Npars);

	CSmoothEmulator emulator(parmap);
	emulator.randy->reset(-time(NULL));
	NPars = emulator.NPars = Npars;

  //cout << "NUMBER OF parameter IS " << NPars <<endl;

	Theta.resize(NPars);
	ThetaTest.resize(ntest);
	for(itest=0;itest<ntest;itest++){
		ThetaTest[itest].resize(NPars);
	}


	// Real function
	real=new CReal_Taylor(NPars,emulator.smooth->MaxRank,emulator.randy);
	real->LAMBDA=emulator.LAMBDA;
	emulator.real=real;
	real->RandomizeA(100.0);

  emulator.SetThetaSimplex();

	CLog::Info("NTrainingPts="+to_string(emulator.NTrainingPts)+"\n");
	emulator.CalcYTrainFromThetaTrain();
	emulator.GenerateASamples();

	for(itest=0;itest<ntest;itest++){
		for(ipar=0;ipar<Npars;ipar++)
			ThetaTest[itest][ipar]=1.0-2.0*emulator.randy->ran();
	}




  for(int ipars = 0; ipars < Npars; ipars++)
  {
    file >> ipars;
    string dirname = "modelruns/run" + std::to_string(ipars);
    string shellcommand = "mkdir -p "+dirname;
    system(shellcommand.c_str());



    for (size_t i = 0; i < NPars; i++) {
      modPar.theta[i] = emulator.ThetaTrain[ipars][i];
    }


    modPar.TranslateTheta_to_x();
      string filename ="modelruns/run" + to_string(ipars)+ "/mod_parameters.txt";
      for (size_t i = 0; i < NPars; i++) {
        fptr = fopen(filename.c_str(),"a");
        fprintf(fptr,"%s %f\n", modPar.priorinfo->parname[i].c_str(),modPar.x[i]);
        fclose(fptr);

    }
      modPar.Print();
      cout << "----------------" <<endl;
  }
    //cout << pInfo <<endl;


  return 0;
}
