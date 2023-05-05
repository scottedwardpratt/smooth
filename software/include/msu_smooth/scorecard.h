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
	unsigned int NTrain,NTest;
	vector<double> ThetaTrain;
	vector<double> ThetaTest;
	vector<double> ProbTest;
	double score,YExp,sigmaYExp;
	void CalcScore(CSmoothEmulator *emulator,CSmooth *smooth,vector<double> &ThetaTestSet,double YExpSet,double sigmaYExpSet);
};

#endif