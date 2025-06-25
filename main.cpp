#include "Matrix_02508971.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <cmath>
#include </home/linux_sys/eigen-3.4.0/Eigen/Dense>

int main(){
    int size = 3;
    double sparsity = 0.0; 

    adv_prog_cw::Matrix_02508971<double> mat(size, size, 0);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> ran(-5, 5);
    Eigen::MatrixXd matrix(size, size);
    std::uniform_real_distribution<double> prob(0, 1); // Probability for sparsity

/*     for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            double a = ran(gen);
            mat(i,j) = a;
            matrix(i,j) = a;
        }
    } */

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                // Ensure diagonal dominance
                int a = 1+ran(gen); // Positive values for stability
                mat(i, j) = a;
                matrix(i,j) = a;
            } else if (prob(gen) > sparsity) { // With (1 - sparsity) probability, make it nonzero
                int b = ran(gen);
                mat(i, j) = b;
                matrix(i,j) = b;
            } else {
                mat(i, j) = 0; // Most off-diagonal elements are zero
                matrix(i,j) = 0;
            }
        }
    }

    //std::cout << matrix << "\n";

    //mat.Out(3);

    // Check 3.1:
/*     auto t1_mul_vec = std::chrono::high_resolution_clock::now();
    mat*=4.2;
    auto t2_mul_vec = std::chrono::high_resolution_clock::now();
    auto time_mul_vec = std::chrono::duration_cast<std::chrono::microseconds>(t2_mul_vec - t1_mul_vec);
    std::cout << "Time for multiplication: " << time_mul_vec.count() << " microseconds \n"; */

    //mat.Out(2);
    
    auto t1_vec = std::chrono::high_resolution_clock::now();
    auto det = mat.Determinant();
    auto t2_vec = std::chrono::high_resolution_clock::now();
    auto time_vec = std::chrono::duration_cast<std::chrono::milliseconds>(t2_vec - t1_vec);
    std::cout << "Time for Det: " << time_vec.count() << " milliseconds\n";

    std::cout << "My det: " << det << "\n";

/*     auto t1_div_vec = std::chrono::high_resolution_clock::now();
    mat/=3.6;
    auto t2_div_vec = std::chrono::high_resolution_clock::now();
    auto time_div_vec = std::chrono::duration_cast<std::chrono::microseconds>(t2_div_vec - t1_div_vec);
    std::cout << "Time for divison: " << time_div_vec.count() << " microseconds \n"; */

    auto t1_vec1 = std::chrono::high_resolution_clock::now();
    adv_prog_cw::Matrix_02508971<double> emp_Mat(50,50,0);
    bool hasDet = mat.Inverse(emp_Mat);
    auto t2_vec1 = std::chrono::high_resolution_clock::now();
    auto time_vec1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2_vec1 - t1_vec1);
    std::cout << "Time for Inv: " << time_vec1.count() << " milliseconds\n"; 

    //emp_Mat.Out(5);

   Eigen::MatrixXd inverse = matrix.inverse();

    double sum_norm = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double diff = mat(i,j) - matrix(i,j);
            sum_norm += std::pow(diff, 2);
        }
    }
    double frobeniusNorm = std::sqrt(sum_norm);

    std::cout << "Difference in inv: " << frobeniusNorm << "\n";
    //std::cout << inverse << "\n";

    double determinant = matrix.determinant();
    std::cout << "Real det: " << determinant << "\n";

    double sum1 = std::pow(determinant,2);
    double sum2 = std::pow(det,2);
    double diff2 = std::abs(sum1-sum2);
    
    double norm = std::sqrt(diff2);
    std::cout << "Difference in det: " << norm << "\n"; 

    return 0;
}