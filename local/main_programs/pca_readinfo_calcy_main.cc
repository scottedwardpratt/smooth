#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "msu_commonutils/parametermap.h"
#include "msu_smooth/master.h"
#include "msu_smooth/PCA.h"

using namespace std;

int main(int argc, char *argv[]) {

  if(argc!=2){
    printf("Usage emulator parameter filename\n");
    exit(1);
  }

  PCA *pca = new PCA(argv[1]);

  pca->ReadPCATransformationInfo();
	
	vector<double> Z,Y;
	
	pca->TransformZtoY(Z,Y,vector<double> &SigmaZ,
	vector<vector<double>> &SigmaY_emulator);

  pca->WriteTransformationInfo();

}
