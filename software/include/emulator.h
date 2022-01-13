#ifndef __EMULATOR_H__
#define __EMULATOR_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "parametermap.h"
#include "misc.h"
#include "randy.h"
#include "constants.h"
#include <list>
#include "train.h"

using namespace std;


class CEmulator{
public:
	CRandy *randy;
	CGSLMatrix_Real *gslmatrix;
	unsigned int NPars,MaxRank,NSmooth,NtrainingPts; // To get variance of Y, you calculate for NSmooth different emulator fits
	vector<vector<double>> trainingthetas;
	vector<double> trainingvalues;
	CSmooth smooth;
	
	CEmulator(unsigned int NPars_set,MaxRank_set);
	void SetNTrainingPts(unsigned int NTrainingPts_set);
	void SetTrainingThetasAndValues(vector<vector<double>> theta_training,vector<double> y_training);
	void TrainEmulator();
	void CalcY_AllSmooth(vector<vector<double>> &trainingthetas,vector<double> &smoothvalues);
	void SaveConfiguration(CSmooth *saved_smooth,CSmooth &smooth,double &value);
	
	void SolveForSmoothPars(unsigned int ismooth);
	
};

#endif