#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/priorinfo.h"
#include "msu_commonutils/log.h"

using namespace std;
int main(int argc,char *argv[]){
	//if(argc!=2){
	//	printf("Usage simplextest priorinfo_filename (to be found inside Info)\n");
	//	exit(1);
	//}
	CparameterMap *parmap=new CparameterMap();
	
	//CPriorInfo *priorinfo=new CPriorInfo("Info/prior_info.txt");
	//printf("check a\n");

	parmap->ReadParsFromFile("parameters/simplex_parameters.txt");
	CSimplexSampler *simplex=new CSimplexSampler(parmap);
	
	simplex->SetThetaType1();
	simplex->WriteModelPars("modelruns");
	return 0;
}
