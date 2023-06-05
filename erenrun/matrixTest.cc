#include <iostream>
#include <Eigen/Dense>
#include <chrono>

using namespace std;
using namespace std::chrono;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {

  for (int i = 1; i < 4; i++) {
    int size = pow(10, i);

    MatrixXd m = MatrixXd::Random(size,size); // sizexsize Matrix filled with random numbers between (-1,1)
    VectorXd v = VectorXd::Random(size);
    std::cout << "For the size:\n" << size << std::endl;


    auto start0 = high_resolution_clock::now();
    VectorXd x = m.colPivHouseholderQr().solve(v);
    auto stop0 = high_resolution_clock::now();

    auto duration0 = duration_cast<microseconds>(stop0 - start0);

    cout << "Time taken by the function colPivHouseholderQr: "<< duration0.count() << " microseconds" << endl;

    auto start1 = high_resolution_clock::now();
    VectorXd y = m.householderQr().solve(v);
    auto stop1 = high_resolution_clock::now();

    auto duration1 = duration_cast<microseconds>(stop1 - start1);

    cout << "Time taken by the function householderQr: "<< duration1.count() << " microseconds" << endl;

    auto start2 = high_resolution_clock::now();
    VectorXd z = m.partialPivLu().solve(v);
    auto stop2 = high_resolution_clock::now();

    auto duration2 = duration_cast<microseconds>(stop2 - start2);

    cout << "Time taken by the function partialPivLu: "<< duration2.count() << " microseconds" << endl;

  }
  return 0;

}
