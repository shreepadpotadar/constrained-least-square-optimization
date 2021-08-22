#ifndef RANDNORMMAT_H
#define RANDNORMMAT_H

#pragma once

#include <random>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>


class RandNormMat{
 public:
  // declaring the class 
  RandNormMat(int r, int c);

  // declaring member functions
  double* GenerateRandomNumbers();
  Eigen::MatrixXd CreateMatrix();

  Eigen::MatrixXd matrix;

 private:
  int num_rows_, num_cols_;

};

#endif