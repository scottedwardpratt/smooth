#ifndef __TRAININGINFO_H__
#define __TRAININGINFO_H__
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

class CSmoothEmulator;
class CSmoothMaster;

class CTrainingInfo{
public:
	CTrainingInfo(int NTrainingPts,CObservableInfo *observableinfo,CPriorInfo *priorinfo);
	int NTrainingPts,NObservables;
	vector<CModelParameters *> trainingpars;
	vector<vector<double>> YTrain,SigmaYTrain;
	vector<CModelParameters *> modelpars;
	void ReadTrainingInfo(string rundirname);
	static CSmoothMaster *smoothmaster;
};

#endif