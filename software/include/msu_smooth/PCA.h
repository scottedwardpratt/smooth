#ifndef __PCA_H__
#define __PCA_H__
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <filesystem>

#include "msu_commonutils/parametermap.h"
#include "msu_smooth/master.h"

class PCA{
public:

  PCA(string parameter_filename);
  int nruns,Nobs;
  Eigen::MatrixXd eigvecs;
	Eigen::VectorXd eigvals;
  vector<vector<double>> Y;
	vector<double> SigmaY,Ybar;
	vector<vector<double>> SigmaY_emulator;
	vector<int> NTrainingList;
	string modelruns_dirname;
	CObservableInfo *observable_info;

  void CalcPCA();
  void WriteZTraining();
  void ReadPCATransformationInfo();
	
	void TranslateZtoY(vector<double> &Z,vector<double> &Y,vector<double> &SigmaZ,vector<vector<double>> &SigmaY_emulator);
};

#endif
