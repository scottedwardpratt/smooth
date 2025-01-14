#ifndef __OBSERVABLE_INFO_H__
#define __OBSERVABLE_INFO_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_smoothutils/misc.h"
#include <vector>
#include <map>

namespace NBandSmooth{

	class CObservableInfo{
	public:
		CObservableInfo(string filename);
		unsigned int NObservables;
		vector<string> observable_name;
		vector<double> ALPHA; // point-by-point uncertainty of model
		map<string,unsigned int> name_map;
		unsigned int GetIPosition(string obsname);  // finds position given name of observable
		string GetName(unsigned int iposition);  // finds name give position
		void ReadObservableInfo(string filename);
		void ReadExperimentalInfo(string filename);
		vector<double> YExp,SigmaExp;
		void PrintInfo();
	};

};

#endif
