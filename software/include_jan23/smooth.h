#ifndef __SMOOTH_H__
#define __SMOOTH_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include <list>

using namespace std;

class CSmooth{
public:
	unsigned int MaxRank,NPars;
	unsigned int NCoefficients;
	vector<vector<unsigned int>> IPar;
	vector<unsigned int> dupfactor;
	vector<unsigned int> rank;
	
	CSmooth();
	CSmooth(unsigned int NPars_set);	
	
	double CalcY(vector<double> &A,vector<double> &Lambda,vector<double> &theta);
	void Copy(CSmooth *smooth);
	double CalcY_Remainder(vector<double> &A,vector<double> &Lambda,vector<double> &theta,unsigned int NTrainingPts);
	static vector<unsigned int> factorial;
};

#endif