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
#include "smooth.h"
#include "gslmatrix.h"

class CEmulator{
public:
	unsigned int NPars,NTrainingPts;
	CRandy *randy;
	CSmooth *smooth;
	CGSLMatrix_Real *gslmatrix;
	vector<vector<double>> M;
	
	CEmulator(unsigned int NPars_set);
	
	void SetSize(vector<double> &A,vector<double> &Lambda);
	void SetA_Zero(vector<double> &A);
	void SetA_RanGauss(double Amag,vector<double> &A); // ~ exp(-A^2/2*Amag^2)
	void SetA_RanSech(double Amag,vector<double> &A); // ~ 1/cosh(A/Amag)
	void SetA_Constant(double Amag,vector<double> &A); // =Amag
	void SetLambda_Constant(double LAMBDA,vector<double> &Lambda);
	void CalcAFromTraining(vector<vector<double>> &thetavec,vector<double> Yvalues,vector<double> &A,vector<double> &Lambda);
	
	void SetNTrainingPts(unsigned int NTrainingPts_set);
	void SetThetaRank1(unsigned int &NTrainingPts,vector<vector<double>> &ThetaTrain);
	void SetThetaRank2(unsigned int &NTrainingPts,vector<vector<double>> &ThetaTrain);
	void TuneA(vector<vector<double>> &thetavec,vector<double> &YTrain,double &Amag,unsigned int nmc,vector<double> &A,vector<double> &Lambda);
	double GetLog_AProb(vector<double> &A,vector<double> Lambda,double Amag);

};

#endif