#pragma once

#include <string>
#include <vector>
#include "Eigen/Sparse"

class LPSolver
{
 public:
  LPSolver(const std::string& text);

  void solve();

  private:
  void parse(const char*);
  void findInitialFeasibleDictionary();

  bool isUnbounded() const;


  Eigen::VectorXd objective_vector;
  Eigen::MatrixXd coefficient_matrix;
};
