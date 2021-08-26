#ifndef LANCZOS_H
#define LANCZOS_H

#pragma once

#include <iostream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>


class Lanczos{
 public:
  //Overload constructor
  Lanczos(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const int max_k);

  // member function that implements Lanczos matrix bidiagonalization algorithm
  Eigen::MatrixXd Solution();

  // matrix B stores the bidiagonalised version of matrix A
  Eigen::MatrixXd mat_b;
  Eigen::MatrixXd mat_a;
  Eigen::VectorXd vec_b;

  int max_k;
};


#endif