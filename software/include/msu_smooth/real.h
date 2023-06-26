#ifndef __REAL_H__
#define __REAL_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>#ifndef __EMULATOR_H__
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
#include <iostream>
#include <Eigen/Dense>
#include "msu_smooth/real.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/master.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

class CSmoothEmulator{
public:
	int iY; // labels observable from observable info
	string observable_name;
	CReal *real;
	Eigen::MatrixXd M;

	double SigmaA0,SigmaAMin,SigmaA,SigmaATrial,MCStepSize,MCSigmaAStepSize,LAMBDA;
	unsigned int NMC;   // NMC is for generating independent samplings of A in TuneA
	unsigned int NASample;
	bool TuneAChooseMCMC,ConstrainA0,CutOffA,UseSigmaYReal,FirstTune;
	vector<vector<double>> ASample;
	vector<double> SigmaASample;
	vector<double> A,ATrial;
	vector<vector<double>> ThetaTrain;
	
	CSmoothEmulator(string observable_name_set);
	
	
	void CalcAFromTraining(vector<double> &AA);
	void PrintA(vector<double> &Aprint);

	void SetThetaTrain();
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

	void GenerateASamples();

	void Init();
	
	static CSmoothMaster *smoothmaster;
	static int NPars;
	static CSmooth *smooth;
	static CparameterMap *parmap;
	static Crandy *randy;
	static unsigned int NTrainingPts;


};

#endif

#include <array>
#include <fstream>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include <list>
#include "msu_smooth/smooth.h"
#include "msu_smooth/simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

using namespace std;

class CReal{
public:
	unsigned int NPars;
	CReal();
	virtual void CalcY(vector<double> &theta,double &Y,double &SigmaY);
	void CalcYTrain(vector<double> &YTrain,vector<double> &SigmaYTrain,int NTrainingPts, vector<vector<double>> ThetaTrain);
};

class CReal_Taylor : public CReal{
public:
	Crandy *randy;
	vector<double> A;
	double LAMBDA;
	CReal_Taylor(unsigned int NPars_Set,int maxrank,Crandy *randy);
	CSmooth *smooth;
	void CalcY(vector<double> &theta,double &Y,double &SigmaY);
	// These are functions for generating fake real models
	void RandomizeA(double SigmaReal);
};


class CReal_EEEK : public CReal{
public:
	Crandy *randy;
	vector<double> A;
	double LAMBDA;
	CReal_EEEK(unsigned int NPars_Set,int maxrank,Crandy *randy);
	CSmooth *smooth;
	void CalcY(vector<double> &theta,double &Y,double &SigmaY);
	//double CalcY_1(vector<double> &A,double LAMBDA,vector<double> &theta);
	// These are functions for generating fake real models
	void RandomizeA(double SigmaReal);
};

#endif
