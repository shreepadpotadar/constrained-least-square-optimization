#include <random>
#include <fstream>
#include <iomanip>
#include <string.h>

#include "lanczos_algo.h"
#include "rand_norm_matrix.h"
#include "lag_mult_estimator.h"


int main(int, char**)
{
    int num_rows;
    int num_cols;
    double sigma;

    std::cout << "Enter the number of rows: " << std::endl;
    std::cin >> num_rows;
    std::cout << "Enter the number of columns: " << std::endl;
    std::cin >> num_cols;
    std::cout << "Enter the value of sigma: " << std::endl;
    std::cin >> sigma;
    
    Eigen::MatrixXd mat_a;
    Eigen::VectorXd vec_b(num_rows);
    
    
    // creates a 2D matrix of size num_rows X num_cols with normal random elements having zero mean and 1/sqrt(num_cols) variance
    RandNormMat randnormmat = RandNormMat(num_rows, num_cols);
    mat_a = randnormmat.CreateMatrix();


    // creates 1D vector of size num_rows X 1 with normal random elements having zero mean and unit variance
    RandNormMat randnormvec = RandNormMat(num_rows, 1);
    double* temp_array = new double[num_rows];
    temp_array = randnormvec.GenerateRandomNumbers();
    for ( int i = 0; i < num_rows; i++) 
    {
        vec_b[i] = sigma * temp_array[i];
    }
    // deallocating the memmory of the temp array
    delete []temp_array;


    // Generates the lower bidigonal matrix of size num_cols X num_cols using Lanczos Bidiagonalization algorithm
    Eigen::MatrixXd mat_b;
    Lanczos lanczos = Lanczos(mat_a, vec_b);
    mat_b = lanczos.Solution();


    // Cholesky decomposition 
    Eigen::MatrixXd L;                                             // stotes the Cholesky decomposition of a matrix
    Eigen::MatrixXd mat_temp;
    Eigen::MatrixXd mat_u(num_rows, num_cols);                     // Matrix U of size (num_rows X num_cols)
    Eigen::MatrixXd mat_identity(num_cols, num_cols);     
    Eigen::MatrixXd mat_u_tilda(num_rows, num_cols);               // Matrix U_tilda of size (num_rows X num_cols)

    mat_identity.setIdentity();                                    // creates Identity matrix of size (num_cols X num_cols)
    mat_temp = Eigen::MatrixXd::Zero(num_rows-num_cols,num_cols);  // defines a zero matrix of size (num_rows-num_cols X num_cols) 

    
    // Cholesky decomposition of (B^T*B-sigma^2 * Identity matrix)
    L = (mat_b.transpose() * mat_b - pow(sigma, 2) * mat_identity).llt().matrixL();
    mat_u.topRows(L.transpose().rows()) = L.transpose();
    mat_u.bottomRows(mat_temp.rows()) = mat_temp;                   // Appends (num_rows-num_cols X num_cols) zero matrix to U matrix 

    
    // creates matrix U_tilda by asigning diagonal elements of Matrix U to zero
    mat_u_tilda = mat_u;
    for (int j = 0; j < num_cols; j++)
    {
        mat_u_tilda(j, j) = 0.0;
    }


    int num_k = 1;
    double eps = 1e-8;
    double lambda_max = 1E2;
    double lambda_step_size = 1E-2;
    double lambda_lk, lambda_uk, lambda_opt;
    
    double mat_atb_norm = (mat_a.transpose() * vec_b).squaredNorm();

    // Computes the lower limit of the Lagrange multiplier lambda of cost function
    LagMultEstimator lag_low = LagMultEstimator(mat_u, sigma, mat_atb_norm);
    lambda_lk = lag_low.FindOptimalLambda(lambda_max, num_k, lambda_step_size, eps);
    std::cout << "Function value at lambda = " << lambda_lk << " is " << lag_low.EstimateFunctionValue(num_k, lambda_lk) << std::endl;


    // Computes the higher limit of the Lagrange multiplier lambda of cost function
    LagMultEstimator lag_high = LagMultEstimator(mat_u_tilda, sigma, mat_atb_norm);
    lambda_uk = lag_high.FindOptimalLambda(lambda_max, num_k, lambda_step_size, eps);
    std::cout << "Function value at lambda = " << lambda_uk << " is " << lag_high.EstimateFunctionValue(num_k, lambda_uk) << std::endl;


    num_k = 2;
    lambda_step_size = 1E-4;
    // Computes the optimium of the Lagrange multiplier lambda of cost function
    LagMultEstimator lag_optim = LagMultEstimator(mat_u, sigma, mat_atb_norm);
    lambda_opt = lag_optim.FindOptimalLambda(lambda_max, num_k, lambda_step_size, eps);
    std::cout << "Function value at lambda = " << lambda_opt << " is " << lag_optim.EstimateFunctionValue(num_k, lambda_opt) << std::endl;

    return 0;
}

