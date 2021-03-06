//#include "stdafx.h"
#include "lag_mult_estimator.h"
#include <iostream>

LagMultEstimator::LagMultEstimator( const Eigen::MatrixXd &U, const double s, const double atb )
{
    // Initialises varibles
    mat_u = U;
    sigma = s;
    mat_atb_norm = atb;

    // checks if num_rows is atleast 1 and num_cols atleast 1
    assert( mat_u.rows() > 0 );
    assert( mat_u.cols() > 0 );
}

double LagMultEstimator::EstimateFunctionValue( const int k, const double lambda )
{
    // check if k is not greater than the num_cols of matrix U
    if (k > mat_u.cols())
    {
        std::cout << "K value is :" << k << std::endl;
    }
    assert( k <= mat_u.cols() );

    Eigen::MatrixXd mat_z(k, k);
    Eigen::MatrixXd identitity_mat(k, k);
    Eigen::VectorXd basis_vec(k);
    
    // initiases to identity matrix
    identitity_mat.setIdentity();
    
    // creates a first basis vector e1 = [1, 0, 0, 0, ...] of dim k by 1
    for ( int i = 0; i < k; i++ )
    {
        if ( i == 0 )
        {
            basis_vec[i] = 1.0;
        }
        else
        {
            basis_vec[i] = 0.0;
        }
    } 

    mat_z = ( mat_u.block(0, 0, k, k).transpose() * mat_u.block(0, 0, k, k) + ( pow( sigma, 2 ) + lambda ) * identitity_mat ).inverse();

    return mat_atb_norm * basis_vec.transpose() * mat_z * mat_z * basis_vec;
}

double LagMultEstimator::FindOptimalLambda(double lambda_max, const double lambda_step, const double eps)
{
    // defines the lambda parameteres
    int k = 1;
    int num_k;
    bool flag = true;
    double lambda_start = - pow(sigma, 2); 
    double lambda_end = lambda_max;
    double target_func_val = mat_u.cols();

    // Find lower value for num_k
    do{
        if (this->EstimateFunctionValue(k, lambda_start) >= target_func_val)
        {
            num_k = k;
            flag = false;
        }
        
        k++;

        if ( k > mat_u.cols())
        {
            num_k = mat_u.cols();
            flag = false;
        }
    }while(flag);

    if (num_k == 1)
    {
        num_k = 2;
    }

    // checks if target is beyond the max lambda range. If so, it increases by 2 fold
    // TO DO: implements this logic inside do-while loop 
    if ( this->EstimateFunctionValue(num_k, lambda_end) > target_func_val )
    {
        lambda_max = 2 * lambda_max;
        lambda_end = lambda_max;
    }


    /* implements the binary search algorithm to find the lambda value
     * corresponds to target function value */
    double lambda_mid;
    bool is_found = false;
    do
    {
        lambda_mid = ( lambda_start + lambda_end ) / 2.0;

        if ( target_func_val < this->EstimateFunctionValue(num_k, lambda_mid))
        {
            lambda_start = lambda_mid;
        }
        else
        {
            lambda_end = lambda_mid;
        }

        if ( ( lambda_end - lambda_start ) <= eps )
        {
            is_found = true;
        }
    } while ( !is_found );
    
    return lambda_start;
}
