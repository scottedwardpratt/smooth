#include "parametermap.h"
#include "constants.h"
#include "smooth.h"
#include "emulator.h"

using namespace std;
int main(int argc,char *argv[]){
	CparameterMap *parmap=new CparameterMap();
	double y,yreal,accuracy,r2,average_accuracy=0.0,sigmay,ybar,y2bar,sigmaybar=0.0;
	long long unsigned int nsigmay=0;
	unsigned int isample,itrain,ntest=1000,ipar,ireal,nreal=10;
	vector<double> Theta;
	double SigmaReal=10.0;

	parmap->ReadParsFromFile("parameters.txt");	
	
	CRandy *randy=new CRandy(-time(NULL));	
	CSmoothEmulator emulator(parmap);
	CReal_Taylor *real= new CReal_Taylor(emulator.NPars, randy);
	emulator.real=real;
	
	Theta.resize(emulator.NPars);
	emulator.randy->reset(-time(NULL));

	emulator.SetThetaSimplex();
	for(ireal=0;ireal<nreal;ireal++){
		accuracy=0.0;
		real->RandomizeA(SigmaReal);
		real->CalcYTrain(emulator.YTrain, emulator.NTrainingPts, emulator.ThetaTrain);

		for(itrain=0;itrain<emulator.NTrainingPts;itrain++){
			yreal=emulator.real->CalcY(emulator.ThetaTrain[itrain]);
		}
		emulator.GenerateASamples();

		for(itrain=0;itrain<emulator.NTrainingPts;itrain++){
			yreal=real->CalcY(emulator.ThetaTrain[itrain]);			
			y=emulator.smooth->CalcY(emulator.ASample[0],emulator.LAMBDA,emulator.ThetaTrain[itrain]);
			if(fabs(yreal-y)>0.001){
				printf("Emulator fails!\n");
				printf("%u, %8.4f, %g, %g\n",itrain,emulator.ThetaTrain[itrain][0],y,yreal);
				exit(1);
			}
		}
		
		for(itrain=0;itrain<ntest;itrain++){
			for(ipar=0;ipar<emulator.NPars;ipar++){
				if(emulator.NPars==1){
					Theta[ipar]=-1.0+(2.0/double(ntest))*(0.5+itrain);
				}
				else{
					do{
						r2=0.0;
						Theta[ipar]=1.0-2.0*emulator.randy->ran();
						r2+=Theta[ipar]*Theta[ipar];
					}while(r2>1.0);
				}
			}
			yreal=real->CalcY(Theta);	
					
			ybar=y2bar=0.0;
			for(isample=0;isample<emulator.NASample;isample++){
				y=emulator.smooth->CalcY(emulator.ASample[isample],emulator.LAMBDA,Theta);
				y2bar+=y*y;
				ybar+=y;
				accuracy+=(y-yreal)*(y-yreal);
			}
			y2bar=y2bar/double(emulator.NASample);
			ybar=ybar/double(emulator.NASample);
			sigmay=sqrt(y2bar-ybar*ybar);
			sigmaybar+=sigmay;
			nsigmay+=1;

		}
		accuracy=sqrt(accuracy/double(emulator.SigmaY0*emulator.SigmaY0*emulator.NASample*ntest));
		accuracy*=100.0;
		average_accuracy+=accuracy;
	}
	sigmaybar=sigmaybar/double(nsigmay);
	average_accuracy=average_accuracy/double(nreal);
//	printf("<accuracy>=%g%%\n",average_accuracy/sqrt(2.0));
	return 0;
}