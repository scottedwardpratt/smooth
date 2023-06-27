#ifndef __SMOOTH_MASTER_H__
#define __SMOOTH_MASTER_H__
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
//#include "msu_smooth/simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>
#include "msu_smooth/real.h"
#include "msu_commonutils/log.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

class CSmoothEmulator;

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
	
	void ReadTrainingInfo(string rundir);
	void GenerateASamples();
	void TuneA();
	void SetThetaTrain();
};

#endif
