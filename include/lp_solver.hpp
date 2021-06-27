#pragma once

#include <string>
#include <vector>
#include "Eigen/Sparse"

class LPSolver
{
 public:
  LPSolver(const std::string& text);

  void Solve();


  private:
  void parse(const char*);
  void translate_to_dictionary_form();

  Eigen::VectorXd objective_vector;
  Eigen::MatrixXd coefficient_matrix;
};
