#ifndef __SMOOTH_MASTER_H__
#define __SMOOTH_MASTER_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/smooth.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

class CSmoothEmulator;
class CTrainingInfo;
class CPriorInfo;
class CObservableInfo;
class CModelParInfo;

class CSmoothMaster{
public:
	int NTrainingPts;
	CSmoothMaster(CparameterMap *parmap_set);
	CparameterMap *parmap;
	unsigned int NPars;
	vector<CSmoothEmulator *> emulator;
	CTrainingInfo *traininginfo;
	CObservableInfo *observableinfo;
	CPriorInfo *priorinfo;
	Crandy *randy;
	CSmooth *smooth;
	string ModelRunDirName;
	
	void ReadTrainingInfo();
	void GenerateASamples();
	void Tune();
	void Tune(string obsname);
	void Tune(int iY);
	void SetThetaTrain();
	void CalcY(int iY,CModelParameters *modelpars,double &Y,double &SigmaY);
	void CalcY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY);
	void CalcAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY);
	void TestAtTrainingPts();
};

#endif
