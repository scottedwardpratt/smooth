#include "parametermap.h"
#include "constants.h"
#include "smooth.h"
#include "emulator.h"
#include "gslmatrix.h"
#include "parameterinfo.h"
#include "math.h"

using namespace std;
int main(int argc,char *argv[]){
	//create real model
	vector<double> ymodel;
	vector<double> yfake, xval, yval;
	for(int i=0;i<1000;i++){
		ymodel[i]=cos(0.1*(double) i/100);
	}
	vector<double> xmin, xmax;
	for(int i=0; i<10;i++){
		xmin[i]=0;
		xmax[i]=10;
	}
	CPriorInfo *info = new CPriorInfo(xmin, xmax);
	
	//generate points and translate them 
	//run simplex, evaluate real y at simplex points and that's training points
	//use that for rest of emulator building
	
	CparameterMap *parmap=new CparameterMap();
	double y,yreal,accuracy,r2,average_accuracy=0.0,sigmay,ybar,y2bar,sigmaybar=0.0;
	long long unsigned int nsigmay=0;
	unsigned int isample,itest,ntest=1000,ipar,ireal,nreal=10;
	vector<double> Theta;
	double SigmaReal = 5.0;

	parmap->ReadParsFromFile("parameters.txt");
//	printf("hello\n");
	CSmoothEmulator emulator(parmap);
//	printf("hello\n");
	
	CSmooth *smooth=new CSmooth(parmap);
//	CSmoothEmulator emulator2(smooth);
	
	Theta.resize(emulator.NPars);
	emulator.randy->reset(-time(NULL));

	emulator.SetThetaSimplex();
	
	
	emulator.YTrain=yfake;
	
	//build real y for emulator
	for(ireal=0;ireal<nreal;ireal++){
		accuracy=0.0;
		//sprintf(filename,"testresults/NPars%d_Lambda%g_TrainType%d_real%u.txt",
		//	emulator.NPars,emulator.LAMBDA,emulator.simplex->TrainType,ireal);
		//fptr=fopen(filename,"w");
		emulator.RandomizeRealA(SigmaReal);
//		emulator.CalcYTrainFromRealA(emulator.YTrain, emulator.NTrainingPts, emulator.ThetaTrain);	
		//training points
		for(int i=0;i<emulator.NTrainingPts;i++){
			yfake[i]=cos(0.1* emulator.CalcYReal(emulator.ThetaTrain[i]));
		}
		//emulator.CalcAFromTraining(emulator.A);
		emulator.GenerateASamples();
		for(itest=0;itest<emulator.NTrainingPts;itest++){
			yreal=emulator.CalcYReal(emulator.ThetaTrain[itest]);
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
			yreal=emulator.CalcYReal(Theta);
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
		accuracy=sqrt(accuracy/double(emulator.SigmaY0*emulator.SigmaY0*emulator.NASample*ntest));
		accuracy*=100.0;
		printf("accuracy=%7.3f%%\n",accuracy);
		average_accuracy+=accuracy;

	}
	
	for(int i=0;i<ntest;i++){
		xval[i]=info->CParameterTranslateTheta_to_x(info->min, info->max, yreal, "range");
		yval[i]=yreal;
	}
	info->Write(xval,yval,"emulator.txt");
	
	return 0;
	
}
