#ifndef __SMOOTH_MCMC_H__
#define __SMOOTH_MCMC_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/master.h"

using namespace std;
using namespace NMSUPratt;
using namespace NBandSmooth;

namespace NBandSmooth{
	
	class CLLCalc;

	class CMCMC{
	public:
		CparameterMap *parmap;
		CSmoothMaster *master;
		Crandy *randy;
		
		CMCMC();
		CMCMC(CSmoothMaster *master);
		unsigned int NPars,NObs,NSkip;
		vector<CModelParameters> trace;
		vector<CModelParameters> burntrace;
		double stepsize;
		void ClearTrace(); // erases trace info so one can start over.
		void ClearBurnTrace();
		void ClearTrace(CModelParameters *modpars); // doesn't start ath theta=0
		void ClearBurnTrace(CModelParameters *modpars);
		
		void BurnInMetropolis(unsigned int Nburn);
		void PerformMetropolisTrace(unsigned int Ntrace);
		void BurnInLangevin(unsigned int Nburn);
		void PerformLangevinTrace(unsigned int Ntrace);
		void Langevin(unsigned int nsteps);
		void WriteTrace(string filename);
		
		void CalcLL(CModelParameters *modpars,double &LL);
		void CalcLLPlusDerivatives(CModelParInfo *modpars,double &LL,vector<double> &dLL_dtheta);
		CLLCalc *llcalc;
	};
	
	class CLLCalc{
	public:
		CLLCalc();
		unsigned int NPars,NObs;
		vector<double> Y,SigmaY,SigmaY_emulator;
		vector<vector<double>> dYdTheta;
		CObservableInfo *obsinfo;
		CPriorInfo *priorinfo;
		CSmoothMaster *master;
		virtual void CalcLL(CModelParameters *modpars,double &LL);
		virtual void CalcLLPlusDerivatives(CModelParInfo *modpars,double &LL,vector<double> &dLL_dtheta);
	};
	
	class CLLCalcSmooth : public CLLCalc{
	public:
		void CalcLL(CModelParameters *modpars,double &LL);
		void CalcLLPlusDerivatives(CModelParInfo *modpars,double &LL,vector<double> &dLL_dtheta);
	};

};

#endif
