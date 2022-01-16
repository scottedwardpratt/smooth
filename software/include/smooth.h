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
#include <list>

using namespace std;
class CSmooth_old;

class CSmooth{
public:
	unsigned int MaxRank,NPars;
	unsigned int NCoefficients;
	vector<vector<unsigned int>> IPar;
	vector<unsigned int> dupfactor;
	vector<unsigned int> rank;
	
	CSmooth();
	CSmooth(unsigned int MaxRank_set,unsigned int NPars_set);
	
	
	double CalcY(vector<double> &A,vector<double> &Lambda,vector<double> &theta);
	void Copy(CSmooth *smooth);
	double CalcY_Remainder(vector<double> &A,vector<double> &Lambda,vector<double> &theta,unsigned int NTrainingPts);
	static vector<unsigned int> factorial;
};

#endif