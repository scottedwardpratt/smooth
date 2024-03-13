#include "msu_smoothutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smoothutils/log.h"
using namespace std;
int main(){
	NMSUUtils::CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile("parameters/emulator_parameters.txt");
	NBandSmooth::CSmoothMaster master(parmap);
	master.ReadCoefficientsAllY();
	NBandSmooth::CModelParameters *modpars=new NBandSmooth::CModelParameters(); // contains info about single point
	modpars->priorinfo=master.priorinfo;
	master.priorinfo->PrintInfo();
	
	//  Calc Observables
	NBandSmooth::CObservableInfo *obsinfo=master.observableinfo;
	vector<double> Y(obsinfo->NObservables);
	vector<double> SigmaY(obsinfo->NObservables);
	master.CalcAllY(modpars,Y,SigmaY);
	cout << "---- EMULATED OBSERVABLES ------\n";
	for(unsigned int iY=0;iY<obsinfo->NObservables;iY++){
		cout << obsinfo->GetName(iY) << " = " << Y[iY] << " +/- " << SigmaY[iY] << endl;
	}

	return 0;
}
