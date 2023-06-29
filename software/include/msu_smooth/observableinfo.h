#ifndef __OBSERVABLE_INFO_H__
#define __OBSERVABLE_INFO_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "msu_commonutils/misc.h"
#include <vector>
#include <map>

class CObservableInfo{
public:
	CObservableInfo(string obs_info_filename_set);
	string observable_info_filename;
	int NObservables;
	vector<string> observable_name,unit;
	vector<double> SigmaA0; // representative spread of coefficients
	map<string,int> name_map;
	int GetIPosition(string obsname);  // finds position given name of observable
	string GetName(int iposition);  // finds name give position
	void ReadObservableInfo(string observable_info_filename);
	void PrintInfo();
};



#endif
