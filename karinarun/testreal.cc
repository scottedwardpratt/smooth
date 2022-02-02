#include "parametermap.h"
#include "constants.h"
#include "smooth.h"
#include "emulator.h"
#include "gslmatrix.h"

using namespace std;

int main(int argc,char *argv[]){
	CReal_Taylor *real_taylor = new CReal_Taylor(3, -time(NULL));
	vector<double> &theta;
	double y;

	theta.resize(real_taylor.NPars);
	for(unsigned int ic=0;ic<real_taylor.RealA.size();ic++){
		real_taylor.RealA[ic]=real_taylor.randy->ran_gauss();
	}
	
	for(int itest=0;itest<real_taylor.NPars;itest++){
		theta[itest]=real_taylor.randy->ran_gauss();
	}
	
	y=real_taylor.CalcY(theta);
	
	
	
	return 0;
}