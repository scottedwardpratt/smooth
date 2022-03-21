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
#include "real.h"

class CSmoothEmulator{
public:
	unsigned int NPars,NTrainingPts;
	CRandy *randy;
	CSmooth *smooth;
//	CReal *real;
	Eigen::MatrixXd M;
	CReal_Taylor *real_taylor;


	double SigmaY0,SigmaYMin,SigmaY,SigmaYTrial,MCStepSize,MCSigmaYStepSize,LAMBDA;
	unsigned int NMC;   // NMC is for generating independent samplings of A in TuneA
	unsigned int NASample;
	bool TuneAChooseMCMC,ConstrainA0,CutOffA;
	vector<vector<double>> ASample;
	vector<double> SigmaYSample;
	vector<double> A,ATrial;
	vector<double> YTrain;
	vector<vector<double>> ThetaTrain;
	CSimplexSampler *simplex;
	
	CSmoothEmulator(CparameterMap *parmap);
	CSmoothEmulator(CSmooth *smooth);
	void CalcYTrainFromThetaTrain();
	void CalcAFromTraining(vector<double> &AA);
	void PrintA(vector<double> &Aprint);
	
	void SetNTrainingPts(unsigned int NTrainingPts_set);
	void SetThetaSimplex();
	void TuneA();
	void TuneAMCMC();
	void TuneAPerfect();
	double GetLog_AProb(vector<double> &AA,double SigmaY);

	void SetA_Zero(vector<double> &A);
	void SetA_RanGauss(double ASigmaY,vector<double> &AA);
	void SetA_Constant(double ASigmaY,vector<double> &AA);
	void SetA_RanSech(double ASigmaY,vector<double> &AA);

	double SigmaYbar;
	int NSigmaY;
	
	CReal *real;
	vector<double> RealA;

	void GenerateASamples();
	
	void Init(CSmooth *smooth);

};

#endif