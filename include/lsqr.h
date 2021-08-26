#ifndef LSQR_H
#define LSQR_H


#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>


class LSQR
{
 public:
  //Overload constructor
  LSQR(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const double eps);

  Eigen::MatrixXd SolveForX();

  Eigen::MatrixXd mat_a;
  Eigen::VectorXd vec_b;
  Eigen::MatrixXd vec_x;

  double eps;

};


#endif