#include "parametermap.h"
#include "constants.h"
#include "smooth.h"
#include "emulator.h"
#include "gslmatrix.h"

using namespace std;
int main(int argc,char *argv[]){
	CparameterMap *parmap=new CparameterMap();
	double y,yreal,accuracy;
	unsigned int isample,itest,ntest=1000,ipar,ireal,nreal=10;
	vector<double> Theta;

	parmap->ReadParsFromFile("parameters.txt");
	CSmoothEmulator emulator(parmap);

	Theta.resize(emulator.NPars);
	emulator.randy->reset(-time(NULL));

	emulator.SetThetaSimplex();
	printf("NTrainingPts=%u\n",emulator.NTrainingPts);

	FILE *fptr;
	char filename[150];
	for(ireal=0;ireal<nreal;ireal++){
		accuracy=0.0;
		sprintf(filename,"testresults_multiD/real%u.txt",ireal);
		fptr=fopen(filename,"w");
		emulator.CalcYTrainFromRealA();	
		//emulator.CalcAFromTraining(emulator.A);
		emulator.GenerateASamples();
		for(itest=0;itest<emulator.NTrainingPts;itest++){
			yreal=emulator.CalcRealYFromRealA(emulator.ThetaTrain[itest]);
			y=emulator.smooth->CalcY(emulator.ASample[1],emulator.Lambda,emulator.ThetaTrain[itest]);
			if(fabs(yreal-y)>0.001){
				printf("Emulator fails!\n");
				printf("%u, %8.4f: %g =? %g\n",itest,emulator.ThetaTrain[itest][0],y,yreal);
			}
		}

		for(itest=0;itest<ntest;itest++){
			for(ipar=0;ipar<emulator.NPars;ipar++){
				Theta[ipar]=1.0-2.0*emulator.randy->ran();
				//Theta[ipar]=-1.0+(2.0/double(ntest))*(0.5+itest);
			}
			yreal=emulator.CalcRealYFromRealA(Theta);
			fprintf(fptr,"%10.7f ",yreal);
			for(isample=0;isample<emulator.NASample;isample++){
				y=emulator.smooth->CalcY(emulator.ASample[isample],emulator.Lambda,Theta);
				accuracy+=(y-yreal)*(y-yreal);
				fprintf(fptr," %10.7f",y);
			}
			fprintf(fptr,"\n");
		}
		fclose(fptr);
		printf("accuracy=%10.3e\n",sqrt(accuracy/double(emulator.NASample*ntest)));
	}
	
	return 0;
}