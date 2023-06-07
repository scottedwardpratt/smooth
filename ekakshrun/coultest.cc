#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/simplex.h"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage smoothy parameter filename (assumed to be found inside parameters/)\n");
		exit(1);
	}
	CparameterMap *parmap=new CparameterMap();
	unsigned int itest,ntest=10,NPars;
	vector<vector<double>> ThetaTest;
	
	string parfilename="parameters/"+string(argv[1]);
	parmap->ReadParsFromFile(parfilename);
	CSmoothEmulator emulator(parmap);
	NPars=emulator.NPars;

	ThetaTest.resize(ntest);
	for(itest=0;itest<ntest;itest++){
		ThetaTest[itest].resize(NPars);
	}
	emulator.randy->reset(-time(NULL));
	
	NAlternativeParameterSampling::GetParsLHC(ntest,NPars,emulator.randy,ThetaTest);
	printf("ThetaTest to start:\n");
	for(int ipar=0;ipar<NPars;ipar++){
		for(int itest=0;itest<ntest;itest++){
			printf("%8.5f ",ThetaTest[itest][ipar]);
		}
		printf("\n");
	}
	NAlternativeParameterSampling::GetParsCoulomb(ntest,NPars,emulator.randy,ThetaTest);
	
	FILE *fptr=fopen("figs/coulomb.txt","w");
	for(int ipar=0;ipar<NPars;ipar++){
		for(itest=0;itest<ntest;itest++){
			fprintf(fptr,"%8.5f ",ThetaTest[itest][ipar]);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
	return 0;
}