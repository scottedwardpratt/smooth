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
#include "smooth.h"
#include "simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

using namespace std;

class CReal{
public:
	CReal();
	virtual double CalcY(vector<double> &theta);
	void CalcYTrain(vector<double> &YTrain,int NTrainingPts, vector<vector<double>> ThetaTrain);
};

class CReal_Taylor : public CReal{
public:
	CRandy *randy;
	unsigned int NPars;
	vector<double> A;
	double LAMBDA;
	CReal_Taylor(unsigned int NPars_Set,CRandy *randy);
	CSmooth *smooth;
	double CalcY(vector<double> &theta);
	// These are functions for generating fake real models
	void RandomizeA(double SigmaReal);
};

#endif