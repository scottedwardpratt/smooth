#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <functional>


using namespace std;

// Function to calculate the force of the harmonic oscillator
double potential_function_3d(const vector<double>& x) {
  double U  = 0.5*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  return U;
}



vector<double> numerical_gradient(const function<double(const vector<double>&) > & f, const vector<double>& point, double eps = 1e-6) {
    std::vector<double> grad(point.size(), 0.0);
    for (size_t i = 0; i < point.size(); ++i) {
        std::vector<double> point_plus_eps = point;
        std::vector<double> point_minus_eps = point;
        point_plus_eps[i] += eps;
        point_minus_eps[i] -= eps;
        grad[i] = (f(point_plus_eps) - f(point_minus_eps)) / (2 * eps);
    }
    return grad;
}

vector<vector<double>> simulate_dynamics(const function<double(const vector<double>&)> & potential_func, double tau, double total_time, const vector<double>& initial_pos) {
    double dt = tau;
    int num_steps = static_cast<int>(total_time / dt);
    vector<vector<double>> traj(num_steps, vector<double>(3, 0.0));
    traj[0] = initial_pos;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(0.0, std::sqrt(2 * dt));

    for (int i = 1; i < num_steps; ++i) {
        vector<double> gradient = numerical_gradient(potential_func, traj[i - 1]);
        vector<double> random_force = {dis(gen), dis(gen), dis(gen)};
        for (size_t j = 0; j < 3; ++j) {
            traj[i][j] = traj[i - 1][j] - gradient[j] * dt + random_force[j];
        }
    }

    return traj;
}

int main() {
    std::vector<double> initial_pos = {2.0, 4.0, 10.0};
    double total_time = 1000.0;
    double tau = 0.001;

    function<double(const vector<double>&)> potential_func = potential_function_3d;
    auto trajectory = simulate_dynamics(potential_func, tau, total_time, initial_pos);

    std::ofstream outFile("trajectory.txt");

    if (outFile.is_open()) {
        for (const auto& pos : trajectory) {
            outFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }


    return 0;
}
