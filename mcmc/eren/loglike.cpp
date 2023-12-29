#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <math.h>

using namespace std;


void GetLogLikelyhood(vector<double> &y,double &y_obs,vector<double> &LL,vector<double> &Der_LL,double sig,vector<double> &d){
  LL.resize(y.size());
  Der_LL.resize(y.size());

  for(int i=0;i<y.size();i++){
		LL[i] = -pow(((y_obs - y[i])/(sqrt(2)*sig)), 2);
    Der_LL[i] = 2*((y_obs - y[i])/(sqrt(2)*sig))*d[i]*1/(sqrt(2)*sig);
  }
}

int main(){
  double sig = 0.07;
  double y_obs = 0.33;
  vector<double> y_emu, y_der, LL, Der_LL;

  std::ofstream outfile("log_values.txt");

  for (double x = -0.99; x < 1; x+= 0.01) {
    y_emu.push_back(1/(M_PI*sqrt(1 - x*x)));
    y_der.push_back(x/(M_PI*pow(1-x*x,3/2)));
  }

  GetLogLikelyhood(y_emu, y_obs, LL, Der_LL, sig, y_der);
  double x = -0.99;
  for (size_t i = 0; i < LL.size(); i++) {

    cout <<x << " " << LL[i] << " "<< Der_LL[i] << endl;

    outfile << x << " " << LL[i]<<" "<< Der_LL[i] << std::endl;
    x += 0.01;
  }
  outfile.close();
  return 0;
}
