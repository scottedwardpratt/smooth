#include "msu_commonutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_commonutils/log.h"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage smoothy emulator parameter filename");
		exit(1);
	}
	Crandy *randy=new Crandy(time(NULL));
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(string(argv[1]));
	
	CSmoothMaster master(parmap);
	master.randy->reset(-time(NULL));
	
	master.ReadCoefficientsAllY();
	
	CModelParameters *modpars=new CModelParameters(master.priorinfo);
	CObservableInfo *obsinfo=new CObservableInfo("Info/Observable_Info.txt");
	int NModelPars=modpars->NModelPars;
	int NObs=obsinfo->NObservables;
	
	//vector<double> theta(NModelPars);
	vector<double> Y(NObs);
	vector<double> SigmaY(NObs);
	for(int ipar=0;ipar<NModelPars;ipar++){
		modpars->X[ipar]=100.0*randy->ran();
		
		printf("%12s=%g\n",master.priorinfo->parname[ipar].c_str(),modpars->X[ipar]);
	}
	modpars->TranslateX_to_Theta();
	
	
	for(int ipar=0;ipar<NModelPars;ipar++){
		printf("Theta[%d]=%g\n",ipar,modpars->Theta[ipar]);
	}
	
	
	master.CalcAllY(modpars,Y,SigmaY);
	for(int iY=0;iY<NObs;iY++){
		printf("%12s=%g\n",(obsinfo->GetName(iY)).c_str(),Y[iY]);
	}

	return 0;
}
