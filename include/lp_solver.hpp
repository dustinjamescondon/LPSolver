#pragma once

#include <string>
#include <vector>
#include "Eigen/Dense"
using namespace Eigen;

class LPSolver
{
public:
  LPSolver(const std::string& text);

  void solve();
  void pivot(size_t entering, size_t leaving);

private:
  MatrixXd A_B() const;
  MatrixXd A_N() const;
  VectorXd x_B() const;
  VectorXd x_N() const;
  VectorXd z_N() const;
  double objective_value() const;

  void parse(const char*);
  void findInitialFeasibleDictionary();

  bool isUnbounded() const;

  size_t num_basic_vars;
  size_t num_non_basic_vars;
  // basis is an n+m sized vector
  Eigen::ArrayXi basis_indices; // contains the indices of the basis variables
  Eigen::ArrayXi non_basis_indices;
  Eigen::VectorXd equational_objective_coefficients;
  Eigen::MatrixXd equational_matrix;
  Eigen::VectorXd b_vector;
};
