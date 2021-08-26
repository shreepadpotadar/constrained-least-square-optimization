#include "rand_norm_matrix.h"


RandNormMat::RandNormMat(int r, int c)
{
    // Initiaze variables
    num_rows = r;
    num_cols = c;

    /* checks if num_rows is atleast 2 and num_cols atleast 1, because
     * generated matrix should satisfy the condition num_rows > num_cols */
    assert( num_rows > 1 );
    assert( num_cols > 0);
}

Eigen::MatrixXd RandNormMat::CreateMatrix()
{
    int num_elements = num_rows * num_rows;
    double* random_nums = new double[num_elements];

    // creates list of normal random variables with zero mean and 1/sqrt(num_cols) variance
    random_nums = this->GenerateRandomNumbers();

    // creates a 2D matrix by mapping 1D array using Map() function
    Eigen::Map<Eigen::MatrixXd> matrix( random_nums, num_rows, num_cols );

    return matrix;
}

double* RandNormMat::GenerateRandomNumbers()
{
    int num_elements = num_rows * num_cols;
    double* random_nums = new double[num_elements];

    for (int i = 0; i < num_elements; i++)
    {
        random_nums[i] = 0.0;
    }

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd; 

    // initialized with seed from previous random device instance
    std::mt19937 gen(rd()); 
    
    //creates a list of normal random numbers with zero mean and 1/sqrt(num_cols) variance
    for(int i = 0; i < num_elements; i++)
    {
        // instance of class std::normal_distribution with specific mean and stddev
        std::normal_distribution<double> d(0, 1/ sqrt(num_cols)); 

        // get random number with normal distribution using gen as random source
        random_nums[i] =  d(gen);
    }
    return random_nums;
}