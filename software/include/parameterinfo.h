#ifndef __PARAMETER_H__
#define __PARAMETER_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "parametermap.h"
#include "misc.h"
#include "randy.h"
#include "constants.h"
#include <list>
#include "smooth.h"
#include "simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>
#include "real.h"


class CParameterInfo{
public:
	string name;
	string filename;
	
	vector<double> min, max;
	double  width, center;
	double theta;
	int param_num;
	
	CParameterInfo(string file);
	CParameterTranslateTheta_to_x(min, max, width, center);
	CParameterTranslateX_to_Theta();
	Print(name, theta);
};

#endif