#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Function to calculate the force of the harmonic oscillator
double potnetial_function_3d(vector<double> &x;) {
  double U  = 0.5*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  return U
}



double numerical_gradient(double f, vector<double> point, double epi = 1e-6){
  vector<double> grad(point.size(), 0.0);
  vector<double> point_plus_epi(point);
  vector<double> point_minus_epi(point);
  for (int i = 0;i< point.size();i ++){
    point_plus_epi[i] =  point[i] + epi;
    point_plus_epi.push_back(point_plus_epi[i]);

    point_minus_epi[i] = point[i] - epi;
    point_minus_epi.push_back(point_minus_epi[i]);

    grad[i] = (potnetial_function_3d(point_plus_epi) -potnetial_function_3d(point_minus_epi))/(2*epi);
    grad.push_back(grad[i]);

  }
  return grad;
}

// Function to run the Langevin dynamics simulation
std::vector<double> simulate_dynamics(double (*potnetial_function_3d)(double, double, double),
double tau, double total_time,
vector<double> initial_pos) {

  int num_steps = static_cast<int>(total_time / tau);

  vector<vector<double>> traj(num_steps, vector<double>(3, 0.0));
  traj[0] = initial_pos;

  // Random number generation setup
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> dis(0.0, std::sqrt(2 * tau));

  for (int i = 1; i < num_steps; ++i) {
    double force = force_func(x[i-1], k,0);
    double random_force = dis(gen);
    x[i] = x[i-1] + (force * tau) + random_force;
  }

  return x;
}



int main() {

  std::vector<double> initial_pos = {1.0,2.0,10.0};
  double total_time = 10000.0;

  std::ofstream outfile("avg_x2_values.txt");

  // Logarithmic spacing of tau
  double tau_min = std::pow(10, -3);
  double tau_max = std::pow(10, 0.01);
  int num_tau_values = 200;
  std::vector<double> taus;
  std::vector<double> avg_x2_values;

  for (int i = 0; i < num_tau_values; ++i) {
    double tau = tau_min * std::pow(tau_max / tau_min, i / (num_tau_values - 1.0));
    taus.push_back(tau);
    std::vector<double> x = simulate_dynamics(harmonic_force, k, tau, total_time, initial_pos);

    double sum_x2 = 0.0;
    for (double xi : x) {
      sum_x2 += xi * xi;
    }
    double avg_x2 = sum_x2 / x.size();
    avg_x2_values.push_back(avg_x2);

    // Write tau and average x^2 to a file
    outfile << tau << " " << avg_x2 << std::endl;

  }

  outfile.close();
  return 0;
}
