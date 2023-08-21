#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdio>

#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/misc.h"
#include <list>
#include "msu_smooth/smooth.h"
#include "msu_smooth/simplex.h"
#include <iostream>
#include <Eigen/Dense>

using namespace std;


class FakeModel {
public:
    int iY;
    string Yname;
    Crandy *randy;
    double coefficient_sin, coefficient_cos, coefficient_exp;
  	vector<double> A;
  	double LAMBDA;
  	CSmooth *smooth;
    unsigned int NPars;

    FakeModel(unsigned int NPars_Set,int maxrank){
    	NPars=NPars_Set;
      randy = new Crandy(time(NULL));
    	smooth = new CSmooth(NPars,maxrank);
      coefficient_sin = 50.0 * randy->ran();
      coefficient_cos = 50.0 * randy->ran();
      coefficient_exp = 50.0 * randy->ran();
    	LAMBDA=2;
    };

    void RandomizeA(double SigmaReal) {
        if (A.size() != smooth->NCoefficients) {
            A.resize(smooth->NCoefficients);
        }
        for (unsigned int ic = 0; ic < A.size(); ic++) {
            A[ic] = SigmaReal * randy->ran_gauss();
        }
    }

    ~FakeModel() {
        delete randy;
    }

    void GetY_1(int iY, string Yname,vector<double> &X, double &Y, double &SigmaY) {
        int ipar, NPars = X.size();
        double Lambda = 2.5, arg = 0.0;

        for (ipar = 0; ipar < NPars; ipar++) {
            arg += randy->ran() * X[ipar] / (2.0 * M_PI);
        }
        Y = coefficient_sin * sin(arg / Lambda) + coefficient_cos * cos(arg / Lambda);
        arg = 0.0;

        for (ipar = 0; ipar < NPars; ipar++) {
            arg += randy->ran() * X[ipar] / (2.0 * M_PI);
        }
        Y += coefficient_exp * exp(arg / (Lambda * NPars));
        SigmaY = 0.0;

        if (iY == 0) {
            printf("-----------------------\n");
        }
        printf("%s  Y=%g +/- %g\n", Yname.c_str(), Y, SigmaY);
    }

    void GetY_2(int iY, string Yname,vector<double> &X,double &Y,double &SigmaY){
    	Y=smooth->CalcY(A,LAMBDA,X);
    	SigmaY=1.0;
      if (iY == 0) {
          printf("-----------------------\n");
      }
      printf("%s  Y=%g +/- %g\n", Yname.c_str(), Y, SigmaY);
    }


};
