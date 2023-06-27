#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"

using namespace std;
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage smoothy parameter filename (assumed to be found inside parameters/)");
		exit(1);
	}
	CparameterMap *parmap=new CparameterMap();
	string parfilename="parameters/"+string(argv[1]);
	parmap->ReadParsFromFile(parfilename);
	CSmoothMaster master(parmap);

	master.randy->reset(-time(NULL));

	printf("Set %d Training Points\n",master.NTrainingPts);

	master.TuneA();
	master.GenerateASamples();


	return 0;
}
