#include <cmath>
#include <iostream>
#include <filesystem>

using namespace std;
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/randy.h"


#include "fakemodcls.cc"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage fakemodel model_parameters_filename\n");
		exit(1);
	}
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(string(argv[1]));



	CObservableInfo *observableinfo=new CObservableInfo("Info/observable_info.txt");
	CPriorInfo *priorinfo=new CPriorInfo("Info/prior_info.txt");
	int NPars,ipar,iY,itrain;
	double Y,SigmaY;
	bool exists;
	vector<double> X;
	FILE *fptr;
	string filename,Yname;
	char parname[300];
	unsigned int NPars_Set;
  int maxrank = 1-NPars_Set;
  double SigmaReal = 10;
	int seed;

	seed = parmap->getI("SmoothEmulator_Seed",time(NULL));


	NPars=priorinfo->NModelPars;
	X.resize(NPars);
	FakeModel fakeModel(NPars_Set, maxrank);
	fakeModel.RandomizeA(SigmaReal);


	exists = true;
	itrain = 0;
	do {

		filename="modelruns/run"+to_string(itrain)+"/mod_parameters.txt";
		exists=filesystem::exists(filename);
		if (exists) {
			fptr = fopen(filename.c_str(), "r");
			for (ipar = 0; ipar < NPars; ipar++) {
				fscanf(fptr, "%s %lf", parname, &X[ipar]);
			}
			fclose(fptr);
			filename = "modelruns/run" + to_string(itrain) + "/obs.txt";
			fptr = fopen(filename.c_str(), "w");
			for (iY = 0; iY < observableinfo->NObservables; iY++) {
				Yname = observableinfo->GetName(iY);
				fakeModel.GetY_1(iY, Yname, X, Y, SigmaY);
				fprintf(fptr, "%s %g %g\n", Yname.c_str(), Y, SigmaY);
			}
			fclose(fptr);
		}
		itrain += 1;
	} while (exists);



	return 0;
}
