#include "parameterinfo.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>

using namespace std;

CPriorInfo::CPriorInfo(string file, vector<double> x_min, vector<double> x_max){

	filename=file;
//	param_num=param;
	FILE *fptr;
	string interval;
	fptr=fopen(filename.c_str(), "r");
	for(int i=0; i<param_num; i++){
		fscanf(fptr, "%s, %s, %lf, %lf", &name[i], &interval, &min[i], &max[i]);
		min[i]=x_min[i];
		max[i]=x_max[i];
		Print(name[i], interval, min[i], max[i]);
	}
//	CParameterTranslateTheta_to_x(min, max, );
	fclose(fptr);
}

CPriorInfo::CPriorInfo(vector<double> x_min, vector<double> x_max){
	min=x_min;
	max=x_max;
	
}

double CPriorInfo::CParameterTranslateX_to_Theta(vector<double> x_min, vector<double> x_max, double x, string interval){
	//for min, max range
	double theta_min=-1.0;
	double theta_max=1.0;
	double center=0;
	double width=1;
	double theta_of_x;
	string str="range";
	
	if(interval==str){
		for(int i=0; i<x_min.size();i++){
			theta_of_x=((x-x_min[i])/(x_max[i]-x_min[i]))*(theta_max-theta_min)+theta_min;
		}
	}
	
	else{
		for(int i=0;i<x_min.size(); i++){
			theta_of_x=exp(-pow(width*(x_max[i]*(x+(center-x_min[i]))),2)/(2*width));
		}
	}
	//for gaussian
	
	return theta_of_x;
}

double CPriorInfo::CParameterTranslateTheta_to_x(vector<double> x_min, vector<double> x_max, double theta, string interval){
	double theta_min=-1.0;
	double theta_max=1.0;
	double center=0;
	double width=1;
	double x_of_theta;
	string str="range";
	
	if(interval==str){
		for(int i=0; i<x_min.size();i++){
			x_of_theta=((theta-theta_min)/(theta_max-theta_min))*(x_max[i]-x_min[i])+x_min[i];
		}
	}
	
	else{
		for(int i=0; i<x_min.size();i++){
			x_of_theta=exp(-pow((width*(theta+(x_min[i]-center))),2)/(2*x_max[i]));
		}
	}
	
	return x_of_theta;
}


void CPriorInfo::Print(string &name, string &type, double &min, double &max){
	for(int i=0;i<param_num;i++){
		printf("Parameter: %s with type %s and range: %lf %lf\n", name.c_str(), type.c_str(), min, max);
		
	}
}

void CPriorInfo::Write(vector<double> &xval,vector<double> &yval, string filename){
	unsigned int ipt,Npts=xval.size();
	double pt;
	FILE *fptr;
	fptr=fopen(filename.c_str(),"w");
	for(ipt=0;ipt<Npts;ipt++){
		fprintf(fptr,"%4.0f %g\n",xval[ipt],yval[ipt]);
	}
	fclose(fptr);
}

CModelParameter::CModelParameter(vector<double> x_val, vector<double> theta_val, CPriorInfo *info_val){
	x=x_val;
	theta=theta_val;
	info=info_val;
}