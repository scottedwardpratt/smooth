#ifndef __PRIORINFO_H__
#define __PRIORINFO_H__

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
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>

class CPriorInfo{
public:
	CPriorInfo(string parinfo_filename);
	int NModelPars;
	string parinfo_filename;
	vector<string> parname,type; // type is gaussian or linear
	vector<double> xmin, xmax;
};



class CModelParameters{
public:
//has actual values of the parameter
	vector<double> X;
	vector<double> Theta;
	CModelParameters(CPriorInfo *priorinfo_set);
	void TranslateTheta_to_X();
	void TranslateX_to_Theta();
	void Print();

	int NModelPars;
	CPriorInfo *priorinfo;
	
};

#endif
