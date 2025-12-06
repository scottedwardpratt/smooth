#ifndef __TPO_H__
#define __TPO_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include <list>
#include <iostream>
#include <Eigen/Dense>
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/modelparinfo.h"
#include <algorithm>
#include <random>

using namespace NMSUUtils;

namespace NBandSmooth{

	class CTPO{
	public:
		bool PLUS1,INCLUDE_LAMBDA_UNCERTAINTY;
		CparameterMap parmap;
		unsigned int NPars,NTrainingPts,NMC;
		string TPO_Method;
		double LAMBDA,ALPHA; //only used for estimating overall uncertainty
		vector<vector<double>> ThetaTrain;
		vector<bool> TrainingPtsFreeze;
		vector<bool> TrainingPtsRead;
		Eigen::MatrixXd I,J,K;
		string FullModelRunsDirName;
		double RSimplex;
		bool FIRSTCALL;
		CPriorInfo *priorinfo;
		Crandy *randy;
		vector<CModelParameters *> modelparameters;
		CTPO();
		void SetThetaTrain(vector<vector<double>> &theta);
		void WriteModelPars();
		
		void CreateTrainingPts();
		void ReadTrainingPts(); // to start from given positions
		void FreezeTrainingPts();
		void SetTrainingPts(); // Sets those points which were not read in from file
			
		void Optimize();
		void Optimize_MC();
		void OptimizeSphere_MC();
		void OptimizeSimplex_MC();
		void SetThetaSimplex(double R);
		void SetThetaSimplexPlus1(double R);
		void SetThetaLatinHyperCube(vector<vector<double>> &theta);
		void SetThetaRandom(vector<vector<double>> &theta);
		
		double GetSigma2Bar(double Lambda,double ALPHA,double &W11);
		void GetC0DDprime(double Lambda,vector<double> &theta1,vector<double> &theta2,double &C0,double &D,double &Dprime);
		void CalcIJK(double Lambda,vector<double> &Rprior);
		void CalcIJK_Gaussian(double Lambda,double beta);
		void GetIiJiKiGaussian(double Rprior,double Lambda,double theta_a,double theta_b,
		double &I,double &Jaterm,double &Jbterm,double &Kabterm);
		void GetIiJiKiUniform(double Rprior,double Lambda,double theta_a,double theta_b,
		double &I,double &Jaterm,double &Jbterm,double &Kabterm);

	};

};


#endif
