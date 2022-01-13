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
#include "randy.h"
#include "constants.h"
#include <list>
//#include "emulator.h"

using namespace std;
class CSmooth_old;

class CSmooth{
public:
	unsigned int MaxRank,NPars;
	unsigned int NCoefficients;
	vector<vector<unsigned int>> IPar;
	vector<double> A;
	vector<unsigned int> dupfactor;
	vector<unsigned int> rank;
	vector<double> Lambda;
	
	CSmooth();
	CSmooth(unsigned int MaxRank_set,unsigned int NPars_set);
	
	void SetA_Zero();
	void SetA_RanGauss(double Amag); // ~ exp(-A^2/2*Amag^2)
	void SetA_RanSech(double Amag); // ~ 1/cosh(A/Amag)
	void SetA_Constant(double Amag); // =Amag
	double CalcY(vector<double> &theta);
	void Copy(CSmooth *smooth);
	void SetLambda(double lambdaset);
	void CalcAFromTraining(vector<vector<double>> &thetatheta,vector<double> Yvalues);
	double CalcY_Remainder(vector<double> &theta,unsigned int NTrainingPts);
	double GetLog_AProb(double Amag);
	static CRandy *randy;
	static vector<unsigned int> factorial;
};

#endif