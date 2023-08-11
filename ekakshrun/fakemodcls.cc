#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdio>
#include "Crandy.h"
using namespace std;
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_smooth/simplex.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/observableinfo.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/randy.h"

using namespace std;


class FakeModel {
private:
    int iY;
    string Yname;
    Crandy *randy;
    double coefficient_sin, coefficient_cos, coefficient_exp;

public:
    FakeModel(int iY, string Yname) : iY(iY), Yname(Yname) {
        randy = new Crandy(iY + time(NULL));
        coefficient_sin = 50.0 * randy->ran();
        coefficient_cos = 50.0 * randy->ran();
        coefficient_exp = 50.0 * randy->ran();
    }

    ~FakeModel() {
        delete randy;
    }

    void GetY(vector<double> &X, double &Y, double &SigmaY) {
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
};
