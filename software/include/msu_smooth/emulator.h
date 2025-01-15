#ifndef __EMULATOR_H__
#define __EMULATOR_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/master.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"
using namespace NMSUUtils;

namespace NBandSmooth{
	class CSmoothMaster;
	class CTrainingInfo;
	class CPriorInfo;
	class CObservableInfo;
	class CModelParInfo;

	class CSmoothEmulator{
	public:
		unsigned int iY; // labels observable from observable info
		string observable_name;

		double LAMBDA,SigmaA,ALPHA;
		double logP;
		bool GPOPTION;
		vector<vector<double>> ThetaTrain,TTrain;
		Eigen::MatrixXd B,Binv;
		Eigen::VectorXd chi;

		CSmoothEmulator(string observable_name_set);
		double GetCorrelation(vector<double> &theta1,vector<double> &theta2);

		void SetThetaTrain();
		void Tune();
		void Tune(double LambdaSet);
		void CalcSigmaA();
		void CalcSigmaALambda();
		void CalcLogP();
		void CalcB();
		
		void GetYAndUncertainty(vector<double> &Theta,double &Y,double &uncertainty);
		
		void WriteCoefficients();
		void ReadCoefficients();
		
		void Init();

		static CSmoothMaster *smoothmaster;
		static unsigned int NPars;
		static CparameterMap *parmap;
		static Crandy *randy;
		static unsigned int NTrainingPts;

	};

};

#endif
