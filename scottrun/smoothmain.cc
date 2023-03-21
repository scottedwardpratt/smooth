#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "smooth.h"
#include "emulator.h"
#include "msu_commonutils/gslmatrix.h"

using namespace std;
int main(int argc,char *argv[]){
	CparameterMap *parmap=new CparameterMap();
	double y,yreal,accuracy,r2,average_accuracy=0.0,sigmay,ybar,y2bar,sigmaybar=0.0;
	long long unsigned int nsigmay=0;
	unsigned int isample,itest,ntest=5,ipar,ireal,nreal=10;
	vector<double> Theta;
	// This plays the role of the "real" model
	CReal_Taylor *real;

	parmap->ReadParsFromFile("parameters.txt");
	CSmoothEmulator emulator(parmap);

	Theta.resize(emulator.NPars);
	emulator.randy->reset(-time(NULL));
	real=new CReal_Taylor(emulator.NPars,emulator.randy);
	real->LAMBDA=emulator.LAMBDA;
	emulator.real=real;

	emulator.SetThetaSimplex();
	printf("Set %d Training Points\n",emulator.NTrainingPts);

	FILE *fptr;
	char filename[150];
	snprintf(filename,sizeof(filename),"testresults/NPars%u_Lambda%g_NTrain%u.txt",
		emulator.NPars,emulator.LAMBDA,emulator.NTrainingPts);
	fptr=fopen(filename,"w");

	for(ireal=0;ireal<nreal;ireal++){
		printf("------ ireal=%d -----\n",ireal);
		accuracy=0.0;
		real->RandomizeA(100.0);
		real->A[0]=0.0;
		for(int ia=0;ia<100;ia++)
			printf("Areal[%d]=%g\n",ia,real->A[0]);
		Misc::Pause();
		// 
		emulator.CalcYTrainFromThetaTrain();
		emulator.GenerateASamples();
		for(itest=0;itest<emulator.NTrainingPts;itest++){
			yreal=emulator.YTrain[itest];
			y=emulator.smooth->CalcY(emulator.ASample[1],emulator.LAMBDA,emulator.ThetaTrain[itest]);
			if(fabs(yreal-y)>0.001){
				printf("Emulator fails!\n");
				printf("%u, %8.4f: %g =? %g\n",itest,emulator.ThetaTrain[itest][0],y,yreal);
			}
		}

		for(itest=0;itest<ntest;itest++){
			for(ipar=0;ipar<emulator.NPars;ipar++){
				if(emulator.NPars==1){
					Theta[ipar]=-1.0+(2.0/double(ntest))*(0.5+itest);
					//fprintf(fptr,"%10.7f ",Theta[ipar]);
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
			fprintf(fptr,"%3u %3u %10.7f ",ireal,itest,yreal);
			ybar=y2bar=0.0;
			for(isample=0;isample<emulator.NASample;isample++){
				y=emulator.smooth->CalcY(emulator.ASample[isample],emulator.LAMBDA,Theta);
				y2bar+=y*y;
				ybar+=y;
				accuracy+=(y-yreal)*(y-yreal);
				//fprintf(fptr,"%10.7f ",y);
			}
			y2bar=y2bar/double(emulator.NASample);
			ybar=ybar/double(emulator.NASample);
			sigmay=sqrt(y2bar-ybar*ybar);
			sigmaybar+=sigmay;
			nsigmay+=1;
			fprintf(fptr,"%10.7f %10.7f\n",ybar,sigmay);
		}
		accuracy=sqrt(accuracy/double(emulator.SigmaY0*emulator.SigmaY0*emulator.NASample*ntest));
		accuracy*=100.0;
		printf("accuracy=%7.3f%%\n",accuracy);
		average_accuracy+=accuracy;
		
	}
	
	sigmaybar=sigmaybar/double(nsigmay);
	printf("<sigma_y>=%g, <SigmaY>=%g\n",sigmaybar,emulator.SigmaYbar/double(emulator.NSigmaY));
	average_accuracy=average_accuracy/double(nreal);
	printf("<accuracy>=%g%%\n",average_accuracy/sqrt(2.0));
	fclose(fptr);
	
	return 0;
}