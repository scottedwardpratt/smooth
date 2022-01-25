#ifndef __SIMPLEX_SAMPLER_H__
#define __SIMPLEX_SAMPLER_H__
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
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

class CSimplexSampler{
public:
	unsigned int NPars,NTrainingPts,TrainType;
	double RTrain;
	CSimplexSampler(CparameterMap *parmap){
		NPars=parmap->getD("Smooth_NPars",0);
		TrainType=parmap->getI("Smooth_TrainType",1);
		RTrain=parmap->getD("Smooth_RTrain",0.9); 
	}
	void SetThetaType1(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaType2(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaType3(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaType4(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);
	void SetThetaSimplex(vector<vector<double>> &ThetaTrain,unsigned int &NTrain);

};

#endif