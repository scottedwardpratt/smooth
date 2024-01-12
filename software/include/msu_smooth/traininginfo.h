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
#include <iostream>
#include <Eigen/Dense>
#include "msu_commonutils/log.h"
#include "msu_smooth/master.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
using namespace NMSUPratt;

namespace NBandSmooth{
	class CSmoothMaster;

	class CTrainingInfo{
	public:
		CObservableInfo *observableinfo;
		CPriorInfo *priorinfo;
		CTrainingInfo(vector<unsigned int> NTrainingList, CObservableInfo *observableinfo,CPriorInfo *priorinfo);
		unsigned int NTrainingPts,NObservables;
		vector<unsigned int> NTrainingList;
		vector<vector<double>> YTrain,SigmaYTrain;
		vector<CModelParameters *> modelpars;
		void ReadTrainingInfo(string rundirname);
		static CSmoothMaster *smoothmaster;
	};

};

#endif
