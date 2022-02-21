#include "parameterinfo.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>

using namespace std;

CParameterInfo::CParameterInfo(string file, int param){

	filename=file;
	param_num=param;
	FILE *fptr;
	fptr=fopen(filename.c_str(), "r");
	for(int i=0; i<param_num; i++){
		fscanf(fptr, "%s, %s, %lf, %lf", &name, %interval, &min[i], &max[i]);

		Print(name, interval, min, max);
	}
	CParameterTranslate(min, max);
	fclose(fptr);
}

CParameterInfo::CParameterTranslateX_to_Theta(vector<double> x_min, vector<double> x_max, double x){
	//for min, max range
	double theta_min=-1.0;
	double theta_max=1.0;
	double center=0;
	double width=1;
	
	for(int i=0; i<x_min.size();i++){
		theta_of_x=((x-x_min[i])/(x_max[i]-x_min[i]))*(theta_max-theta_min)+theta_min;
	}
	
	//for gaussian
	for(int i=0;i<x_min.size(); i++){
		theta_of_x=exp(-(x+(center-x_min))^2/(2*width));
	}
}

CParameterInfo::CParameterTranslateTheta_to_x(vector<double> x_min, vector<double> x_max, double theta){
	double theta_min=-1.0;
	double theta_max=1.0;
	double center=0;
	double width=1;
	
	//for interval
	for(int i=0; i<x_min.size();i++){
		x_of_theta=((theta-theta_min[i])/(theta_max[i]-theta_min[i]))*(x_max[i]-x_min[i])+x_min[i];
	}
	//for gaussian
	for(int i=0; i<x_min.size();i++){
		theta_of_x=exp(-(x+(x_min-center))^2/(2*x_max));
	}
	
}


CParameterInfo::Print(string &name, string &type, double &min, double &max){
	printf("Parameter: %s with type %s and range: %lf %lf\n", name, type, min, max);
}
