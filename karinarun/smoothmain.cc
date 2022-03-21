#include "parametermap.h"
#include "constants.h"
#include "smooth.h"
#include "emulator.h"
#include "gslmatrix.h"

using namespace std;
int main(int argc,char *argv[]){
	CparameterMap *parmap=new CparameterMap();
	double y,yreal,accuracy,r2,average_accuracy=0.0,sigmay,ybar,y2bar,sigmaybar=0.0;
	long long unsigned int nsigmay=0;
	unsigned int isample,itest,ntest=1000,ipar,ireal,nreal=10;
	vector<double> Theta;
	double SigmaReal = 5.0;

	parmap->ReadParsFromFile("parameters.txt");	
	
	CRandy *randy=new CRandy(-time(NULL));	
	CSmoothEmulator emulator(parmap);
	//	CReal_Taylor *real_taylor= new CReal_Taylor(emulator.NPars, randy);
	
	CSmooth *smooth=new CSmooth(parmap);
	//	CSmoothEmulator emulator2(smooth);
	
	Theta.resize(emulator.NPars);
	emulator.randy->reset(-time(NULL));

	emulator.SetThetaSimplex();
//	cout << "NTrainingPts is:" << emulator.NTrainingPts << endl;
	
	//FILE *fptr;
	//char filename[150];
	for(ireal=0;ireal<nreal;ireal++){
		accuracy=0.0;
		//sprintf(filename,"testresults/NPars%d_Lambda%g_TrainType%d_real%u.txt",
		//	emulator.NPars,emulator.LAMBDA,emulator.simplex->TrainType,ireal);
		//fptr=fopen(filename,"w");
		emulator.real_taylor->RandomizeRealA(SigmaReal);
		emulator.real_taylor->CalcYTrainFromRealA(emulator.YTrain, emulator.NTrainingPts, emulator.ThetaTrain);	
//			cout << "YTrain is:" << emulator.Yrain.size() << endl;
		//emulator.CalcAFromTraining(emulator.A);
		emulator.GenerateASamples();
		for(itest=0;itest<emulator.NTrainingPts;itest++){
			yreal=emulator.real_taylor->CalcYReal(emulator.ThetaTrain[itest]);			
			y=emulator.smooth->CalcY(emulator.ASample[0],emulator.LAMBDA,emulator.ThetaTrain[itest]);
			cout << "thetatrain is:" << emulator.ThetaTrain[itest][0] << endl;
//			cout << "y is:" << y << endl;
			if(fabs(yreal-y)>0.001){
				printf("Emulator fails!\n");
				printf("%u, %8.4f, %g, %g\n",itest,emulator.ThetaTrain[itest][0],y,yreal);
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
			yreal=emulator.real_taylor->CalcYReal(Theta);	
					
			//fprintf(fptr,"%10.7f ",yreal);
			ybar=y2bar=0.0;
			for(isample=0;isample<emulator.NASample;isample++){
				y=emulator.smooth->CalcY(emulator.ASample[isample],emulator.LAMBDA,Theta);
				y2bar+=y*y;
				ybar+=y;
				accuracy+=(y-yreal)*(y-yreal);
				//fprintf(fptr," %10.7f",y);
			}
			y2bar=y2bar/double(emulator.NASample);
			ybar=ybar/double(emulator.NASample);
			sigmay=sqrt(y2bar-ybar*ybar);
			sigmaybar+=sigmay;
			nsigmay+=1;

			//fprintf(fptr,"\n");
		}
		//fclose(fptr);
		cout << accuracy << endl;
		accuracy=sqrt(accuracy/double(emulator.SigmaY0*emulator.SigmaY0*emulator.NASample*ntest));
		accuracy*=100.0;
	//	printf("accuracy=%7.3f%%\n",accuracy);
		average_accuracy+=accuracy;
	}
	sigmaybar=sigmaybar/double(nsigmay);
//	printf("<sigma_y>=%g, <SigmaY>printf=%g\n",sigmaybar,emulator.SigmaYbar/double(emulator.NSigmaY));
	average_accuracy=average_accuracy/double(nreal);
//	printf("<accuracy>=%g%%\n",average_accuracy/sqrt(2.0));
	return 0;
}