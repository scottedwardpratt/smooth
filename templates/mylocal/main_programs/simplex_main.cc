#include "msu_commonutils/parametermap.h"
#include "msu_smooth/simplex.h"
#include "msu_commonutils/log.h"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage simplextest simplex_parameters_filename\n");
		exit(1);
	}
	NMSUPratt::CparameterMap *parmap=new CparameterMap();

	parmap->ReadParsFromFile(string(argv[1]));
	NBandSmooth::CSimplexSampler *simplex=new NBandSmooth::CSimplexSampler(parmap);
	simplex->SetThetaSimplex();
	simplex->WriteModelPars();
	return 0;
}
