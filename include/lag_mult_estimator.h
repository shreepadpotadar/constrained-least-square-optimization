#ifndef LAGMULTESTIMATOR_H
#define LAGMULTESTIMATOR_H

#pragma once

#include <math.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>


class LagMultEstimator{
 public:
  //Overload constructor
  LagMultEstimator( const Eigen::MatrixXd &U, const double s, const double atb );

  double EstimateFunctionValue( const int k, const double lambda );

  double FindOptimalLambda( double max_lambda, const double lambda_step, const double eps );

  Eigen::MatrixXd mat_u;
  double sigma;
  double mat_atb_norm;

};

#endif