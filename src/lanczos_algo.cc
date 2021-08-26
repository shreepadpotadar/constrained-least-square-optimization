#include "lanczos_algo.h"


Lanczos::Lanczos(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const int max_k)
{
    /***************************************
     * Initialize the variables
     ***************************************/
    mat_a = A;
    vec_b = b;
    this->max_k = max_k;

    /* checks if num_rows is atleast 2 and num_cols atleast 1, because
     * generated matrix should satisfy the condition num_rows > num_cols */
    assert( mat_a.rows() > 2 );
    assert( mat_a.cols() > 1 );
}

Eigen::MatrixXd Lanczos::Solution()
{
    //Lanczos bidiagonalization algorithm
    int num_rows = mat_a.rows();
    int num_cols = mat_a.cols();
    //int min_dim = fmin(num_rows, num_cols);
    int min_dim = max_k;
    int num_of_elements = (min_dim + 1) * min_dim;

    // initialises the intermediate vectors used in the algorithm
    Eigen::VectorXd delta(min_dim);
    Eigen::VectorXd gamma(min_dim);
    Eigen::VectorXd vec_p(num_rows);
    Eigen::VectorXd vec_q(num_cols);
    vec_q.setZero();
    Eigen::VectorXd vec_w_even(num_rows);
    Eigen::VectorXd vec_w_odd(num_cols);
    

    int k;    // counter that keeps track of the column number
    int i = 1;
    bool flag = true;
    vec_w_even = vec_b;
    do{
        k = (i + 1 ) / 2;
        // checks if i is even or not
        if (i % 2 != 0)
        {
            // executes if i is odd number
            if (i==1)
            {
                vec_p = vec_w_even / vec_w_even.norm();
                vec_w_odd = mat_a.transpose() * vec_p - vec_w_even.norm() * vec_q;
            }
            else
            {
                vec_p = vec_w_even / delta[k-2];
                vec_w_odd = mat_a.transpose() * vec_p - delta[k-2] * vec_q; 
            }
            gamma[k-1] = vec_w_odd.norm();
        }
        else
        {
            // executes if i is even number
            vec_q = vec_w_odd / gamma[k-1];
            vec_w_even = mat_a * vec_q - gamma[k-1] * vec_p;
            delta[k-1] = vec_w_even.norm();
        } 
        
        if (i > 2)
        {
            // checks if the stipping criteria is true
            if ( gamma[k-1] * delta[k-2] == 0 )
            {
                flag = false;
            }
            else if( k == min_dim and i % 2 == 0)
            {
                flag = false;
            }
        }
    
        i++;

    }while(flag);

    double* data = new double[(int) num_of_elements];

    /* creates 1D array that can be transformed into a 2D lower bidiagonal matrix 
     * from gamma and delta vectors */
    k = 0;    // counter keeps track of the location of element stored in 1-D array
    for (int j = 0; j < min_dim; j++)
    {
        for (int i = 0; i <= min_dim; i++)
        {
            if (i == j)
            {
                data[k] = gamma[i];
            }
            else if (i == j + 1)
            {
                data[k] = delta[j];
            }
            else
            {
                data[k] = 0.0;
            }
            k++;
        }
    }

    // creates a 2D matrix by mapping 1D array using Map() functio
    Eigen::Map<Eigen::MatrixXd> mat_b(data, (min_dim + 1), min_dim);

    return mat_b;
}