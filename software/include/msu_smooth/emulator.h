#ifndef __EMULATOR_H__
#define __EMULATOR_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include <list>
#include "msu_smooth/smooth.h"
#include "msu_smooth/simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>
#include "msu_smooth/real.h"
#include "msu_commonutils/log.h"

class CSmoothEmulator{
public:
	unsigned int NPars,NTrainingPts;
	Crandy *randy;
	CSmooth *smooth;
	CReal *real;
	Eigen::MatrixXd M;

	double SigmaA0,SigmaAMin,SigmaA,SigmaATrial,MCStepSize,MCSigmaAStepSize,LAMBDA;
	unsigned int NMC;   // NMC is for generating independent samplings of A in TuneA
	unsigned int NASample;
	bool TuneAChooseMCMC,ConstrainA0,CutOffA,UseSigmaYReal,FirstTune;
	vector<vector<double>> ASample;
	vector<double> SigmaASample;
	vector<double> A,ATrial;
	vector<double> YTrain,SigmaYTrain;
	vector<vector<double>> ThetaTrain;
	CSimplexSampler *simplex;

	CSmoothEmulator(CparameterMap *parmap);
	CSmoothEmulator(CSmooth *smooth);
	void CalcYTrainFromThetaTrain();
	void CalcYTrainFromThetaTrain_EEEK();
	void CalcAFromTraining(vector<double> &AA);
	void PrintA(vector<double> &Aprint);

	void InitTrainingPtsArrays(unsigned int NTrainingPts_set);
	void SetThetaSimplex();
	void TuneA();
	void TuneAMCMC();
	void TuneAMCMC_withSigma();
	void TuneAPerfect();
	double GetLog_AProb(vector<double> &AA,double SigmaA);

	void SetA_Zero(vector<double> &A);
	void SetA_RanGauss(double ASigmaA,vector<double> &AA);
	void SetA_Constant(double ASigmaA,vector<double> &AA);
	void SetA_RanSech(double ASigmaA,vector<double> &AA);

	double SigmaAbar;
	int NSigmaA;

	vector<double> RealA;

	void GenerateASamples();

	void Init(CparameterMap *parmap);

};

#endif
