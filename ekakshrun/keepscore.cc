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
	double YExp,SigmaYExp,SigmaYReal,average_score_yexp=0.0,average_score_real=0.0;
	unsigned int itest,ntest=100,ipar,NPars,ireal,nreal=10;
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

	emulator.SetThetaSimplex();
	CLog::Info("NTrainingPts="+to_string(emulator.NTrainingPts)+"\n");

	for(itest=0;itest<ntest;itest++)

	{
		for(ipar=0;ipar<NPars;ipar++)
		{
			ThetaTest[itest][ipar]=1.0-2.0*emulator.randy->ran();
		}

	}

	real=new CReal_Taylor(NPars,emulator.smooth->MaxRank,emulator.randy);
	real->LAMBDA=emulator.LAMBDA;
	emulator.real=real;
	for(ireal=0;ireal<nreal;ireal++)
	{

		printf("------ ireal=%d -----\n",ireal);
		real->RandomizeA(100.0);

		emulator.CalcYTrainFromThetaTrain();
		emulator.GenerateASamples();

		// Choose a value for YExp -- just some value of YReal for a random theta.
		SigmaYReal=0.0;
		SigmaYExp=0.05;
		average_score_yexp=0.0;
		for(int iyexp=0;iyexp<20;iyexp++)
		{
			for(ipar=0;ipar<NPars;ipar++)
			{
				Theta[ipar]=-1.0+2.0*emulator.randy->ran();
			}
			real->CalcY(Theta,YExp,SigmaYReal);
			scorecard.CalcScore(&emulator,ThetaTest,YExp,SigmaYExp);
			printf("score=%g\n",scorecard.score);
			average_score_yexp += scorecard.score;

		}
		average_score_yexp = average_score_yexp/20;
		printf("average_score_yexp=%g\n",average_score_yexp);

		average_score_real += average_score_yexp/double(nreal);


	}
	printf("average_score = %g\n", average_score_real);
	return 0;

}
