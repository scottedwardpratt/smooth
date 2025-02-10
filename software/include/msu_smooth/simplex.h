#ifndef __SIMPLEX_SAMPLER_H__
#define __SIMPLEX_SAMPLER_H__
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

	class CSimplexSampler{
	public:
		CparameterMap parmap;
		unsigned int NPars,NTrainingPts,TrainType;
		double LAMBDA; //only used for estimating overall uncertainty
		vector<vector<double>> ThetaTrain;
		Eigen::MatrixXd I,J,K;
		string ModelDirName;
		double RGauss,RGauss1,RGauss2;
		CPriorInfo *priorinfo;
		Crandy *randy;
		vector<CModelParameters *> modelparameters;
		CSimplexSampler();
		void SetThetaTrain(vector<vector<double>> &theta);
		void Optimize();
		void Optimize_MC();
		void OptimizeSphere_MC();
		void OptimizeSimplex();
		void OptimizeSimplexPlusCenter();
		
		double GetSigma2Bar(double LAMBDA,double ALPHA,double &detB,double &W11);
		void GetC0DDprime(double LAMBDA,vector<double> &theta1,vector<double> &theta2,double &C0,double &D,double &Dprime);
		void CalcIJK(double LAMBDA,double beta);

	};


	namespace NAlternativeParameterSampling{
		// Latin Hyer Cube parameters
		void GetParsLHC(unsigned int NRuns,unsigned int NPars,Crandy *randy,vector<vector<double>> &Theta);
		void GetParsLHC_Modified(unsigned int NRuns,unsigned int NPars,Crandy *randy,vector<vector<double>> &Theta);
		double GetPEShuffle(vector<vector<double>> x);
		// These are used for Coulomb force generated parameters
		void GetParsCoulomb(unsigned int NRuns,unsigned int NPars,Crandy *randy,vector<vector<double>> &Theta);
		void CalcEnergy(vector<vector<double>> &x,vector<vector<double>> &vv,vector<vector<double>> &v,double &PE,double &KE,double &Etot);
		void Propagate(double dt,vector<vector<double>> &x,vector<vector<double>> &xx,vector<vector<double>> &v,vector<vector<double>> &vv,double &PE);
		void GetForcePotential(vector<double> &x,vector<double> &xx,vector<double> &F,double &potential);
		//
		void GetParsCoulombHO(Crandy *randy,vector<vector<double>> &Theta);
		void GetFRelVRelHO(vector<double> &x,vector<double> &xx,vector<double> &Frel,double &Vrel);
		void PropagateHO(double dt,vector<vector<double>> &x,vector<vector<double>> &xx,vector<vector<double>> &v,vector<vector<double>> &vv,double &PE);
		void CalcEnergyHO(vector<vector<double>> &x,vector<vector<double>> &vv,vector<vector<double>> &v,double &PE,double &KE,double &Etot);
		void GetFextVextHO(vector<double> &x,vector<double> &Fext,double &Vext);
	};

};


#endif