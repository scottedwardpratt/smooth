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
#include "simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

class CSmoothEmulator{
public:
	unsigned int NPars,NTrainingPts;
	CRandy *randy;
	CSmooth *smooth;
	Eigen::MatrixXd M;

	double Amag,AmagTrial,MCstepsize,LAMBDA,RTrain;
	vector<double> Lambda;
	unsigned int NMC;   // NMC is for generating independent samplings of A in TuneA
	unsigned int NASample,TrainRank;
	vector<vector<double>> ASample;
	vector<double> AmagSample;
	vector<double> A,ATrial;
	vector<double> YTrain;
	vector<vector<double>> ThetaTrain;
	CSimplexSampler *simplex;
	
	CSmoothEmulator(CparameterMap *parmap);
	void CalcAFromTraining(vector<double> &AA);
	void PrintA(vector<double> &Aprint);
	
	void SetNTrainingPts(unsigned int NTrainingPts_set);
	void SetThetaSimplex();
	void TuneA();
	double GetLog_AProb(vector<double> &AA,double Amag);

	void SetA_Zero(vector<double> &A);
	void SetA_RanGauss(double AAmag,vector<double> &AA);
	void SetA_Constant(double AAmag,vector<double> &AA);
	void SetA_RanSech(double AAmag,vector<double> &AA);

	void SetLambda_Constant(double LAMBDA_set);

	// These are functions for generating fake real models
	void CalcYTrainFromRealA();
	double CalcRealYFromRealA(vector<double> &theta); // This also sets up a random RealA
	vector<double> RealA;

	void GenerateASamples();

};

#endif