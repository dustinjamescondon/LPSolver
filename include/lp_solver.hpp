#pragma once

#include <string>
#include <vector>
#include "Eigen/Dense"
using namespace Eigen;

class LPSolver
{
public:
  LPSolver(const std::string& text);

  double solve();
  void pivot(size_t entering, size_t leaving);

private:
  MatrixXd A_B() const;
  MatrixXd A_N() const;
  VectorXd calcX_B() const;
  VectorXd calcZ_N() const;
  VectorXd c_B() const;
  VectorXd c_N() const;
  VectorXd deltaX(size_t entering_index);

  double calcHighestIncrease(unsigned int entering_index, unsigned int& leaving_index_out);
  double objective_value();

  size_t chooseEnteringVariable() const;
  void findInitialFeasibleDictionary();

  bool isPrimalFeasible() const;
  bool isDualFeasible() const;
  bool isUnbounded();
  bool isOptimal();

  size_t num_basic_vars;
  size_t num_non_basic_vars;
  // basis is an n+m sized vector
  Eigen::ArrayXi basis_indices; // contains the indices of the basis variables
  Eigen::ArrayXi non_basis_indices;
  Eigen::MatrixXd equational_matrix;
  Eigen::VectorXd c_vector;
  Eigen::VectorXd z_vector;
  Eigen::VectorXd b_vector;
  Eigen::VectorXd x_vector;
};
