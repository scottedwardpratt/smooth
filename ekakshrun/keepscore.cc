#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/smooth.h"
#include "msu_smooth/emulator.h"
#include "msu_commonutils/gslmatrix.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/scorecard.h"
#include "msu_smooth/master.h"
#include "msu_smooth/traininginfo.h"

using namespace std;
int main(int argc,char *argv[]){
	CScoreCard scorecard;
	if(argc!=2){
		printf("Usage smoothy parameter filename (assumed to be found inside parameters/)\n");
		exit(1);
	}
	double average_score_yexp=0.0,average_score_real=0.0;
	unsigned int ntest=100,ipar,NPars, itest;
	vector<vector<double>> ThetaTest;
	vector<double> Theta;


	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(string(argv[1]));
	CSmoothMaster master(parmap);
	CSimplexSampler *simplex=new CSimplexSampler(parmap);
	master.randy->reset(-time(NULL));
	NPars=master.NPars;


	Theta.resize(NPars);
	ThetaTest.resize(ntest);

	master.ReadTrainingInfo();

	for(itest=0;itest<ntest;itest++){
		ThetaTest[itest].resize(NPars);
	}

	simplex->SetThetaSimplex();

	for(itest=0;itest<ntest;itest++)
	{
		for(ipar=0;ipar<NPars;ipar++)
		{
			ThetaTest[itest][ipar]=1.0-2.0*master.randy->ran();
		}
	}

	for(itest = 0; itest < ntest; itest++)
	{
		for(ipar=0;ipar<NPars;ipar++)
		{
			Theta[ipar]=-1.0+2.0*master.randy->ran();
		}
		double YExp = master.traininginfo->YTrain[0][itest];
		double SigmaYExp = master.traininginfo->SigmaYTrain[0][itest]; // we can replace 0 with the correct index for the observable
		scorecard.CalcScore(master.emulator[0],ThetaTest,YExp,SigmaYExp);
		printf("score=%g\n",scorecard.score);
		average_score_yexp += scorecard.score;
	}
	return 0;
}
