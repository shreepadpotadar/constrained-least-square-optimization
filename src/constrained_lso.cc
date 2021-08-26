#include <random>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>

#include <time.h>

#include "lanczos_algo.h"
#include "rand_norm_matrix.h"
#include "lag_mult_estimator.h"
#include "lsqr.h"



int main(int, char**)
{
    /******************************************************************************
     * Takes input from the user through the console
     ******************************************************************************/
    int num_rows;
    int num_cols;
    double sigma;
    std::cout << "Enter the number of rows: " << std::endl;
    std::cin >> num_rows;
    std::cout << "Enter the number of columns: " << std::endl;
    std::cin >> num_cols;
    std::cout << "Enter the value of sigma: " << std::endl;
    std::cin >> sigma;
    std::cout << '\n';



    
    /******************************************************************************
     * Generage matrix A of size MxN (num_rows, num_cols) and vector b of size  Mx1 
     * Matrix A: A(i, j) has zero mean and variance of 1/ sqrt(num_cols)
     * Vector b: b[i] = sigma * X_i where X_i is unit normal variable
     ******************************************************************************/
    // Start measuring time
    clock_t time_0 = clock();

    Eigen::MatrixXd mat_a;
    Eigen::VectorXd vec_b(num_rows);

    RandNormMat randnormmat = RandNormMat(num_rows, num_cols);
    mat_a = randnormmat.CreateMatrix();
    
    RandNormMat randnormvec = RandNormMat(num_rows, 1);
    double* temp_array = new double[num_rows];
    temp_array = randnormvec.GenerateRandomNumbers();
    for ( int i = 0; i < num_rows; i++) 
    {
        vec_b[i] = sigma * temp_array[i];
    }
    // deallocating the memmory of the temp array
    delete []temp_array;

    // Stop measuring time and calculate the elapsed time
    clock_t time_1 = clock();
    double elapsed = double(time_1 - time_0)/CLOCKS_PER_SEC;
    std::cout << "Time taken for data creation: " << elapsed << " seconds.\n" << std::endl;


    /******************************************************************************
     * Lanczos Bidiagonalisation Algorithm
     ******************************************************************************/
    Eigen::MatrixXd mat_b;
    int max_k = 5;
    Lanczos lanczos = Lanczos(mat_a, vec_b, max_k);
    mat_b = lanczos.Solution();



    /*******************************************************************************
     * Cholesky Decomposition
     ******************************************************************************/
    Eigen::MatrixXd L;                                             // stotes the Cholesky decomposition of a matrix
    Eigen::MatrixXd mat_temp;
    Eigen::MatrixXd mat_u(num_rows, max_k);                     // Matrix U of size (num_rows X max_k)
    Eigen::MatrixXd mat_identity(max_k, max_k);     
    Eigen::MatrixXd mat_u_tilda(num_rows, max_k);               // Matrix U_tilda of size (num_rows X max_k)
    mat_identity.setIdentity();                                    // creates Identity matrix of size (max_k X max_k)
    mat_temp = Eigen::MatrixXd::Zero(num_rows-max_k, max_k);  // defines a zero matrix of size (num_rows-max_k X max_k) 
    
    // Cholesky decomposition of (B^T*B-sigma^2 * Identity matrix)
    L = (mat_b.transpose() * mat_b - pow(sigma, 2) * mat_identity).llt().matrixL();
    mat_u.topRows(L.transpose().rows()) = L.transpose();
    mat_u.bottomRows(mat_temp.rows()) = mat_temp;                   // Appends (num_rows-max_k X max_k) zero matrix to U matrix 
    
    // creates matrix U_tilda by asigning diagonal elements of Matrix U to zero
    mat_u_tilda = mat_u;
    for (int j = 0; j < max_k; j++)
    {
        mat_u_tilda(j, j) = 0.0;
    }

    // Stop measuring time and calculate the elapsed time
    clock_t time_2 = clock();
    elapsed = double(time_2 - time_1)/CLOCKS_PER_SEC;
    std::cout << "Time taken for Lanczos algo and cholesky decomposition: " << elapsed << " seconds.\n" << std::endl; 



    /******************************************************************************
     * Computes the optimum Lagrange multiplier of a cost function
     ******************************************************************************/
    int num_k = 1;
    double eps = 1e-8;
    double lambda_max = 1E4;
    double lambda_step_size = 1E-2;
    double lambda_lk, lambda_uk, lambda_opt;
    double mat_atb_norm = (mat_a.transpose() * vec_b).squaredNorm();

    // Computes the lower limit of the Lagrange multiplier lambda of cost function
    LagMultEstimator lag_low = LagMultEstimator(mat_u, sigma, mat_atb_norm);
    lambda_lk = lag_low.FindOptimalLambda(lambda_max, num_k, lambda_step_size, eps);
    std::cout << "Lower limit of the Lagrange multiplier: " << lambda_lk << std::endl;
    std::cout << '\n';

    // Computes the higher limit of the Lagrange multiplier lambda of cost function
    LagMultEstimator lag_high = LagMultEstimator(mat_u_tilda, sigma, mat_atb_norm);
    lambda_uk = lag_high.FindOptimalLambda(lambda_max, num_k, lambda_step_size, eps);
    std::cout << "Upper limit of the Lagrange multiplier: " << lambda_uk << std::endl;
    std::cout << '\n';

    // Computes the optimium of the Lagrange multiplier lambda of cost function
    num_k = 2;
    lambda_step_size = 1E-4;
    LagMultEstimator lag_optim = LagMultEstimator(mat_u, sigma, mat_atb_norm);
    lambda_opt = lag_optim.FindOptimalLambda(lambda_max, num_k, lambda_step_size, eps);
    std::cout << "Optimum value of the Lagrange multiplier: " << lambda_opt << std::endl;
    std::cout << '\n';


    // Stop measuring time and calculate the elapsed time
    clock_t time_3 = clock();
    elapsed = double(time_3 - time_2)/CLOCKS_PER_SEC;
    std::cout << "Time taken for finding optimal lambda: " << elapsed << " seconds.\n" << std::endl; 

    /******************************************************************************
     * LSQR: Least Square QR Factorisation Algorithm
    ******************************************************************************/
    Eigen::MatrixXd mat_a_tilda(num_rows + num_cols, num_cols);
    Eigen::MatrixXd mat_i;
    mat_i.setIdentity(num_cols, num_cols);
    mat_a_tilda.topRows(mat_a.rows()) = mat_a;
    mat_a_tilda.bottomRows(mat_i.rows()) = sqrt(lambda_opt) * mat_i; 

    Eigen::VectorXd vec_b_tilda(num_rows + num_cols, 1);
    for (int i = 0; i < vec_b_tilda.size(); i++)
    {
        if ( i < num_rows)
        {
            vec_b_tilda[i] = vec_b[i];
        }
        else
        {
            vec_b_tilda[i] = 0;
        }
        
    }
    
    // Estiamtion of vec_x using LSQR algorithm
    Eigen::MatrixXd vec_x;
    LSQR lsqr = LSQR(mat_a, vec_b, eps);
    vec_x = lsqr.SolveForX();



    /*******************************************************************************
     * Comparision of results of numerical and analytical solution of the cost
     * funciton of a simplified random optimization
     ******************************************************************************/
    std::cout << "Error estimated by numerical solution: " << std::endl;
    std::cout << (0.5*(mat_a * vec_x - vec_b).squaredNorm() / num_cols) << std::endl;
    std::cout << '\n';

    std::cout << "Error estimated by analytical solution: " << std::endl;
    std::cout << 0.5 * pow( (sqrt( ((1.0 * num_rows)/num_cols) * ( 1 + pow(sigma, 2) ) ) - 1 ), 2) << std::endl;
    std::cout << '\n'; 

    // Stop measuring time and calculate the elapsed time
    clock_t time_4 = clock();
    elapsed = double(time_4 - time_3)/CLOCKS_PER_SEC;
    std::cout << "Time taken for computing X: " << elapsed << " seconds.\n" << std::endl; 


    // Total time taken to execute the program
    elapsed = double(time_4 - time_0)/CLOCKS_PER_SEC;
    std::cout << "Total time taken: " << elapsed << " seconds.\n" << std::endl;  

    return 0;
}

