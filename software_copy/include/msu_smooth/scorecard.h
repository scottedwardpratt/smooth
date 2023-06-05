#ifndef __SCORECARD_H__
#define __SCORECARD_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include <list>
#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

class CScoreCard{
public:
	double score,YExp,SigmaYExp;
	void CalcScore(CSmoothEmulator *emulator,vector<vector<double>> &ThetaTest,double YExpSet,double SigmaYExpSet);
};

#endif