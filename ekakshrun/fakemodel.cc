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

#include "fakemodels/fake1.cc"

using namespace std;
int main(){
	CObservableInfo *observableinfo=new CObservableInfo("Info/observable_info.txt");
	CPriorInfo *priorinfo=new CPriorInfo("Info/prior_info.txt");
	int NPars,ipar,iY,itrain;
	double Y,SigmaY;
	bool exists;
	vector<double> X;
	FILE *fptr;
	string filename,Yname;
	char parname[300];



	NPars=priorinfo->NModelPars;
	X.resize(NPars);




	vector<FakeModel> fakeModels;
	for (int i = 0; i < observableinfo->NObservables; i++) {
		string Yname = observableinfo->GetName(i);
		fakeModels.emplace_back(i, Yname);
	}

	exists = true;
	itrain = 0;
	do {
		// ... other code ...

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
				fakeModels[iY].GetY(X, Y, SigmaY); // Use the FakeModel's GetY method
				fprintf(fptr, "%s %g %g\n", Yname.c_str(), Y, SigmaY);
			}
			fclose(fptr);
		}
		itrain += 1;
	} while (exists);

	return 0;
}
