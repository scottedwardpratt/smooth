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
#include "train.h"

using namespace std;

class CSmooth{
public:
	static CRandy *randy;
	double Amag,A0,LAMBDA;
	vector<double> A1;
	vector<vector<double>> A2;
	vector<vector<vector<double>>> A3;
	vector<vector<vector<vector<double>>>> A4;
	vector<vector<vector<vector<vector<double>>>>> A5;
	vector<int> dupfactor1;
	vector<vector<int>> dupfactor2;
	vector<vector<vector<int>>> dupfactor3;
	vector<vector<vector<vector<int>>>> dupfactor4;
	vector<vector<vector<vector<vector<int>>>>> dupfactor5;
	int MaxRank=5,NPars;
	vector<long long int> factorial;
	CSmooth(unsigned int MaxRank_set,unsigned int NPars_set,unsigned int NTrainingPoints_set);
	void SetLambda(double lambdaset){
		LAMBDA=lambdaset;
	}
	void SetA_Zero();
	void SetA_RanGauss(double Amag); // ~ exp(-A^2/2*Amag^2)
	void SetA_RanSech(double Amag); // ~ 1/cosh(A/Amag)
	void SetA_Constant(double Amag); // =Amag
	double CalcY(vector<double> &theta);
	CSmoothEmulator *smoothemulator;
};

class CSmooth_Linear{
public:
	static CRandy *randy;
	unsigned int NCoefficients;
	unsigned int NTrainingPoints;
	unsigned int NPars;
	vector<vector<unsigned int>> IPar;
	vector<double> A;
	vector<unsigned int> dupfactor;
	vector<unsigned int> rank;
	CSmooth_Linear(unsigned int MaxRank_set,unsigned int NPars_set,unsigned int NTrainingPoints_set);
	void Import(CSmooth *smooth);
	void Export(CSmooth *smooth);
	void void SetA_Zero();
	void SetA_RanGauss(double Amag); // ~ exp(-A^2/2*Amag^2)
	void SetA_RanSech(double Amag); // ~ 1/cosh(A/Amag)
	void SetA_Constant(double Amag); // =Amag
	
	
};
	

#endif