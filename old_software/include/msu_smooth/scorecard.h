#ifndef __SCORECARD_H__
#define __SCORECARD_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/smooth.h"

class CScoreCard{
public:
	double score,YExp,SigmaYExp,SigmaYReal;
	int itest,isample,NTest;
	double yi,Pi,Pibar,Pi2bar;
	void CalcScore(CSmoothEmulator *emulator,vector<vector<double>> &ThetaTest,double YExpSet,double SigmaYExpSet);
};

#endif
