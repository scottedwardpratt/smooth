#ifndef __SMOOTH_H__
#define __SMOOTH_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smoothutils/log.h"
#include <list>

using namespace std;
using namespace NMSUPratt;

namespace NBandSmooth{

	class CSmooth{
	public:
		unsigned int MaxRank,NPars;
		unsigned int NCoefficients;
		vector<unsigned int> factorial;
		vector<vector<unsigned int>> IPar;
		vector<unsigned int> dupfactor;
		vector<unsigned int> rank;
		bool UseRFactor;

		CSmooth();
		CSmooth(unsigned int NPars_Set,unsigned int maxrank);
		CSmooth(CparameterMap *parmap);
		void InitArrays();
		void Copy(CSmooth *smooth);

		double CalcY(vector<double> &A,double LAMBDA,vector<double> &theta);
		double CalcY_FromMtot(vector<double> &A,vector<double> &Mtot);
		double CalcY_Remainder(vector<double> &A,double LAMBDA,vector<double> &theta,unsigned int NTrainingPts);
		double CalcY_Remainder_FromMtot(vector<double> &A,unsigned int NTrainingPts,vector<double> &Mtot);
		double GetRFactor(double LAMBDA,vector<double> &theta);
		double GetT(unsigned int ic,double LAMBDA,vector<double> &theta);
		
		void CalcYDYDTheta(vector<double> &A,double lambda,vector<double> &theta,double &Y,vector<double> &dYdTheta);

	};

};

#endif
