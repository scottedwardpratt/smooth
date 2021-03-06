#ifndef __SMOOTH_H__
#define __SMOOTH_H__
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

class CSmooth{
public:
	unsigned int MaxRank,NPars;
	unsigned int NCoefficients;
	vector<unsigned int> factorial;
	vector<vector<unsigned int>> IPar;
	vector<unsigned int> dupfactor;
	vector<unsigned int> rank;
	bool UseRFactor;

	
	CSmooth();
	CSmooth(unsigned int NPars_Set);
	CSmooth(CparameterMap *parmap);
	void InitArrays();
	
	double CalcY(vector<double> &A,double LAMBDA,vector<double> &theta);
	void Copy(CSmooth *smooth);
	double CalcY_Remainder(vector<double> &A,double LAMBDA,vector<double> &theta,unsigned int NTrainingPts);
	double GetRFactor(double LAMBDA,vector<double> &theta);
	double GetM(int ic,double LAMBDA,vector<double> &theta);
	
};

#endif