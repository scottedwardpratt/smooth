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
  int nruns;
  Eigen::MatrixXd eigvals, eigvecs;
  vector<vector<double>> Y,SigmaY;
	vector<int> NTrainingList;
	string modelruns_dirname;

  void CalcPCA();
  void WriteZTraining();
  void ReadPCA();
};

#endif
