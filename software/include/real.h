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

using namespace std;

class CReal{
public:
	CReal();
	virtual double CalcY(vector<double> &theta);
};

class CReal_Taylor : public CReal{
public:
	CRandy *randy;
	unsigned int NPars;
	vector<double> RealA;
	double LAMBDA;
	CReal_Taylor(unsigned int NPars_Set,CRandy *randy);
	CSmooth *smooth;
	double CalcY(vector<double> &theta);
	// These are functions for generating fake real models
	void RandomizeRealA(double SigmaReal);
};

#endif