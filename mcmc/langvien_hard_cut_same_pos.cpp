#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <functional>

using namespace std;

// Function to calculate the potential of the system
double potential_function_nd(const vector<double>& x) {
    if (x[0] <= -1 || x[0] >= 1) {
        return INFINITY;  // Return infinity for hard cutoffs
    } else {
        double U = 0.5 * (x[0] * x[0]);  // 1D potential function
        return U;
    }
}

vector<double> numerical_gradient(const function<double(const vector<double>&)>& f, const vector<double>& point, double eps = 1e-6) {
    vector<double> grad(point.size(), 0.0);
    for (size_t i = 0; i < point.size(); ++i) {
        vector<double> point_plus_eps = point, point_minus_eps = point;
        point_plus_eps[i] += eps;
        point_minus_eps[i] -= eps;
        grad[i] = (f(point_plus_eps) - f(point_minus_eps)) / (2 * eps);
    }
    return grad;
}

vector<vector<double>> simulate_dynamics(const function<double(const vector<double>&)>& potential_func, double tau, double total_time, const vector<double>& initial_pos) {
    double dt = tau;
    int num_steps = static_cast<int>(total_time / dt);
    vector<vector<double>> traj(num_steps, vector<double>(1, 0.0));
    traj[0] = initial_pos;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dis(0.0, sqrt(2 * dt));

    for (int i = 1; i < num_steps; ++i) {
        vector<double> gradient = numerical_gradient(potential_func, traj[i - 1]);
        double random_force = dis(gen);
        double new_x = traj[i - 1][0] - gradient[0] * dt + random_force;

        // Check for hard cutoffs and update position accordingly
        if (new_x > -1 && new_x < 1) {
            traj[i][0] = new_x;
        } else {
            traj[i][0] = traj[i - 1][0];  // Stay at the current position if outside bounds
        }
    }

    return traj;
}

int main() {
    vector<double> initial_pos = {0.5};  // Start within bounds
    double total_time = 10000.0;
    double tau = 0.01;

    function<double(const vector<double>&)> potential_func = potential_function_nd;
    auto trajectory = simulate_dynamics(potential_func, tau, total_time, initial_pos);

    ofstream outFile("trajectory_hard_cut.txt");
    if (outFile.is_open()) {
        for (const auto& pos : trajectory) {
            outFile << pos[0] << endl;
        }
        outFile.close();
    } else {
        cerr << "Unable to open file for writing." << endl;
    }

    return 0;
}
