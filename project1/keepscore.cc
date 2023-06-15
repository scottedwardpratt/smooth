#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/scorecard.h"

using namespace std;
int main(int argc,char *argv[]){
	CScoreCard scorecard;
	if(argc!=2){
		printf("Usage smoothy parameter filename (assumed to be found inside parameters/)\n");
		exit(1);
	}
	CparameterMap *parmap=new CparameterMap();
	double YExp,SigmaYExp,SigmaYReal;
	unsigned int itest,ntest=1000,ipar,NPars;
	vector<vector<double>> ThetaTest;
	vector<double> Theta;

	// This plays the role of the "real" model
	CReal_Taylor *real;

	string parfilename="parameters/"+string(argv[1]);
	parmap->ReadParsFromFile(parfilename);
	CSmoothEmulator emulator(parmap);
	emulator.randy->reset(-time(NULL));
	NPars=emulator.NPars;

	Theta.resize(NPars);
	ThetaTest.resize(ntest);
	for(itest=0;itest<ntest;itest++){
		ThetaTest[itest].resize(NPars);
	}

	// Real function
	real=new CReal_Taylor(NPars,emulator.smooth->MaxRank,emulator.randy);
	real->LAMBDA=emulator.LAMBDA;
	emulator.real=real;
	real->RandomizeA(100.0);


	emulator.SetThetaSimplex();
	CLog::Info("NTrainingPts="+to_string(emulator.NTrainingPts)+"\n");
	emulator.CalcYTrainFromThetaTrain();
	emulator.GenerateASamples();

	for(itest=0;itest<ntest;itest++){
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTest[itest][ipar]=1.0-2.0*emulator.randy->ran();
	}

	// Choose a value for YExp -- just some value of YReal for a random theta.
	SigmaYExp=0.25;
	SigmaYReal=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		Theta[ipar]=-1.0+2.0*emulator.randy->ran();
	real->CalcY(Theta,YExp,SigmaYReal);

	scorecard.CalcScore(&emulator,ThetaTest,YExp,SigmaYExp);
	printf("score=%g\n",scorecard.score);

	cout << "pointer val:" << parmap << endl;
	cout << "Npars" << parmap->NPars << endl;


	//cout << "*pointer val:" << parmap << endl;

	return 0;
}
