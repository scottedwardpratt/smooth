#include "parametermap.h"
#include "constants.h"
#include "smooth.h"

using namespace std;
int main(int argc,char *argv[]){
	double F,Fbar=0.0;
	vector<double> theta;
	int maxrank=5,NPars=10,ntries=1000000,itry;
	CRandy *randy=new CRandy(-time(NULL));
	printf("Enter NPars: ");
	scanf("%d",&NPars);
	theta.resize(NPars);
	CSmooth *smooth=new CSmooth(maxrank,NPars,randy);
	smooth->SetLambda(3.0);
	//smooth->SetA_Constant(2.0);
	smooth->SetA_RanGauss(1.0);
	
	for(itry=0;itry<ntries;itry++){
		for(int ipar=0;ipar<NPars;ipar++){
			//theta[ipar]=1.0;
			theta[ipar]=randy->ran_gauss();
		}
		F=smooth->CalcY(theta);
		Fbar+=F;
		//printf("F=%g\n",F);
	}
	Fbar=F/double(ntries);
	printf("<F>=%g\n",F);
	return 0;
}