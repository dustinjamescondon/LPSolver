#pragma once

#include <string>
#include <vector>
#include "Eigen/Dense"
using namespace Eigen;

class LPSolver
{
public:
  LPSolver(const char* filename);
  LPSolver() {}

  void solve();
  void pivot(size_t entering, size_t leaving);

private:
  struct HighestIncreaseResult {
    bool unbounded;
    size_t index;
    double maxIncrease;
  };

  struct LPResult {
    double optimal_val;
    bool isUnbounded;
    bool isInfeasible;
  };

  MatrixXd A_B() const;
  MatrixXd A_N() const;
  VectorXd calcX_B() const;
  VectorXd calcZ_N() const;
  VectorXd c_B() const;
  VectorXd c_N() const;
  VectorXd deltaX(size_t entering_index);

  double primalObjectiveValue() const;
  double dualObjectiveValue() const;

  HighestIncreaseResult calcHighestIncrease(unsigned int entering_index);
  HighestIncreaseResult calcHighestIncrease_Dual(unsigned int leaving_index);

  size_t chooseEnteringVariable() const;
  size_t chooseLeavingVariable_Dual() const;
  void findInitialFeasibleDictionary();

  bool isPrimalFeasible() const;
  bool isDualFeasible();
  bool isUnbounded();
  bool isOptimal();

  void primalSolve(Eigen::VectorXd const& obj_coeff_vector);
  LPResult dualSolve(Eigen::VectorXd const& obj_coeff_vector);

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
