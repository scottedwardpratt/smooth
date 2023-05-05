#ifndef __SIMPLEX_SAMPLER_H__
#define __SIMPLEX_SAMPLER_H__
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
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

class CSimplexSampler{
public:
	unsigned int NPars,NTrainingPts,TrainType;
	double RTrain;
	CSimplexSampler(CparameterMap *parmap){
		NPars=parmap->getD("SmoothEmulator_NPars",0);
		TrainType=parmap->getI("Simplex_TrainType",1);
		RTrain=parmap->getD("Simplex_RTrain",0.9); 
	}
	void SetThetaType1(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaType2(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaType3(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaType4(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaSimplex(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);

};

namespace NAlternativeParameterSampling{
	// Latin Hyer Cube parameters
	void GetParsLHC(int NRuns,int NPars,Crandy *randy,vector<vector<double>> &Theta);
	void GetParsLHC_Modified(int NRuns,int NPars,Crandy *randy,vector<vector<double>> &Theta);
	double GetPEShuffle(vector<vector<double>> x);
	// These are used for Coulomb force generated parameters
	void GetParsCoulomb(int NRuns,int NPars,Crandy *randy,vector<vector<double>> &Theta);
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

#endif