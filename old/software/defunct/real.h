#ifndef __REAL_H__
#define __REAL_H__

#include <array>
#include <fstream>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include <list>
#include "msu_smooth/smooth.h"
#include "msu_smooth/simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace NMSUPratt;

namespace NBandSmooth{

	class CReal{
	public:
		unsigned int NPars;
		CReal();
		virtual void CalcY(vector<double> &theta,double &Y,double &SigmaY);
		void CalcYTrain(vector<double> &YTrain,vector<double> &SigmaYTrain,unsigned int NTrainingPts, vector<vector<double>> ThetaTrain);
	};

	class CReal_Taylor : public CReal{
	public:
		Crandy *randy;
		vector<double> A;
		double LAMBDA;
		CReal_Taylor(unsigned int NPars_Set,unsigned int maxrank,Crandy *randy);
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

};

#endif
