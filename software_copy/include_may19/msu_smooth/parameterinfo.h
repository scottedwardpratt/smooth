#ifndef __PRIOR_H__
#define __PRIOR_H__
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
#include "smooth.h"
#include "simplex.h"
//#include "gslmatrix.h"
#include <iostream>
#include <Eigen/Dense>
#include "real.h"


class CPriorInfo{
public:
	vector<string> name;
	string filename;
	
	vector<double> min, max;
	double  width, center;
	double theta;
	int param_num;
	
	CPriorInfo(string file, vector<double> x_min, vector<double> x_max);
	CPriorInfo(vector<double> x_min, vector<double> x_max);
	double CParameterTranslateTheta_to_x(vector<double> min, vector<double> max,double theta, string interval);
	double CParameterTranslateX_to_Theta(vector<double> x_min, vector<double> x_max, double x, string interval);
	void Print(string &name, string &type, double &min, double &max);
	void Write(vector<double> &xval, vector<double> &yval,string filename);
//	ReadInfo();
};



class CModelParameter{
public:
//has actual values of the parameter
//each object tells you the point if you run 100 points will have 100 objects
	CModelParameter(vector<double> x_val, vector<double> theta_val, CPriorInfo *info_val);
	vector<double> x;
	vector<double> theta;
	CPriorInfo *info;




};
#endif
