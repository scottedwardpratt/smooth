#ifndef __REAL_H__
#define __REAL_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "parametermap.h"
#include "misc.h"
#include "constants.h"
#include "randy.h"
#include <list>


using namespace std;

class CReal{
public:
	CReal();
	virtual void CalcRealY(vector<double> &y, vector<double> &theta);

};

/*class CSmoothEmulator : public CReal{
	void CalcRealY(vector<double> &y, vector<double> &theta);
	double CalcRealYFromRealA(vector<double> &theta);
	void CalcYTrainFromRealA();
	void RandomizeRealA();
};
*/
class CReal_Taylor : public CReal{
	unsigned int NPars;
	CReal_Taylor(unsigned int NPars_Set);
	void CalcRealY(vector<double> &y, vector<double> &theta);
};

#endif