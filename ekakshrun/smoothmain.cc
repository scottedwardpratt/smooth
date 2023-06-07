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
	double y,yreal,accuracy,r2,average_accuracy=0.0,average_expected_accuracy=0.0,sigmay2,ybar,y2bar;
	unsigned int isample,itest,ntest=25,ipar,ireal,nreal=10;
	vector<double> Theta;
	// This plays the role of the "real" model
	CReal_Taylor *real;

	string parfilename="parameters/"+string(argv[1]);
	printf("parfilename=%s\n",parfilename.c_str());
	parmap->ReadParsFromFile(parfilename);
	printf("check aa\n");
	CSmoothEmulator emulator(parmap);
		printf("check a: MaxRank=%d\n",emulator.smooth->MaxRank);

	Theta.resize(emulator.NPars);
	emulator.randy->reset(-time(NULL));
	printf("check b: MaxRank=%d\n",emulator.smooth->MaxRank);
	real=new CReal_Taylor(emulator.NPars,emulator.smooth->MaxRank,emulator.randy);
	real->LAMBDA=emulator.LAMBDA;
	emulator.real=real;

	printf("check c : MaxRank=%d\n",emulator.smooth->MaxRank);

	emulator.SetThetaSimplex();
	printf("Set %d Training Points\n",emulator.NTrainingPts);

	//FILE *fptr;
	//char filename[150];
	//snprintf(filename,150,"testresults/NPars%u_Lambda%g_NTrain%u.txt",
	//	emulator.NPars,emulator.LAMBDA,emulator.NTrainingPts);
	//fptr=fopen(filename,"w");

	for(ireal=0;ireal<nreal;ireal++){
		printf("------ ireal=%d -----\n",ireal);
		accuracy=0.0;
		real->RandomizeA(100.0);
		real->A[0]=0.0;
		//
		emulator.CalcYTrainFromThetaTrain();
		emulator.GenerateASamples();
		for(itest=0;itest<emulator.NTrainingPts;itest++)
		{
			yreal=emulator.YTrain[itest];
			y=emulator.smooth->CalcY(emulator.ASample[1],emulator.LAMBDA,emulator.ThetaTrain[itest]);
			if(fabs(yreal-y)>0.001)
			{
				printf("Emulator fails!\n");
				printf("%u, %8.4f: %g =? %g\n",itest,emulator.ThetaTrain[itest][0],y,yreal);
			}
		}

		accuracy=sigmay2=0.0;
		for(itest=0;itest<ntest;itest++)
		{
			for(ipar=0;ipar<emulator.NPars;ipar++)
			{
				if(emulator.NPars==1)
				{
					Theta[ipar]=-1.0+(2.0/double(ntest))*(0.5+itest);
				}
				else
				{
					do
					{
						r2=0.0;
						Theta[ipar]=1.0-2.0*emulator.randy->ran();
						r2+=Theta[ipar]*Theta[ipar];
					}
					while(r2>1.0);
				}
			}
			yreal=real->CalcY(Theta);
			ybar=y2bar=0.0;
			for(isample=0;isample<emulator.NASample;isample++){
				y=emulator.smooth->CalcY(emulator.ASample[isample],emulator.LAMBDA,Theta);
				y2bar+=y*y;
				ybar+=y;
			}
			y2bar=y2bar/double(emulator.NASample);
			ybar=ybar/double(emulator.NASample);
			printf("ybar=%g =? %g\n",ybar,yreal);

			accuracy+=(ybar-yreal)*(ybar-yreal);
			sigmay2+=y2bar-ybar*ybar;
		}
		accuracy=accuracy/ntest;
		sigmay2=sigmay2/ntest;
		accuracy=sqrt(accuracy);
		sigmay2=sqrt(sigmay2/2.0);
		printf("accuracy=%7.3f, expected accuracy=%7.3f\n",accuracy,sigmay2);
		average_accuracy+=accuracy;
		average_expected_accuracy+=sigmay2;

	}
	average_accuracy=average_accuracy/double(nreal);
	average_expected_accuracy=average_expected_accuracy/double(nreal);

	printf("<accuracy>=%g, <expected accuracy>=%g\n",average_accuracy,average_expected_accuracy);

	return 0;
}
