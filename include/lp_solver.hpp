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
  enum  class State { Optimal, Unbounded, Infeasible };

  struct HighestIncreaseResult {
    bool unbounded;
    size_t index;
    double maxIncrease;
  };

  struct LPResult {
    double optimal_val;
    State state;
  };

  MatrixXd A_B() const;
  MatrixXd A_N() const;
  VectorXd calcX_B() const;
  VectorXd calcZ_N() const;
  VectorXd c_B() const;
  VectorXd c_N() const;
  VectorXd deltaX(size_t entering_index);
  VectorXd deltaZ(size_t leaving_index);

  VectorXd solveA_B_transpose_x_equals_b(VectorXd const& b) const;
  VectorXd solveA_B_x_equals_b(VectorXd const& b) const;

  double primalObjectiveValue() const;
  double dualObjectiveValue(VectorXd const& obj_coeff_vector) const;
  void printOptimalVariableAssignment() const;

  HighestIncreaseResult calcHighestIncrease(unsigned int entering_index);
  HighestIncreaseResult calcHighestIncrease_Dual(VectorXd const& delta_z);

  size_t chooseEnteringVariable_blandsRule() const;
  size_t chooseEnteringVariable_bland() const;
  size_t chooseDualLeavingVariable_largestCoeff() const;
  size_t chooseDualLeavingVariable_blandsRule() const;
  void findInitialFeasibleDictionary();

  bool isPrimalFeasible() const;
  bool isDualFeasible();
  bool isOptimal();

  LPResult primalSolve();
  LPResult dualSolve(Eigen::VectorXd const& obj_coeff_vector);

  void updateSubmatrices() const;

  /*--------------------------------------------------
   * mutable because these are basically used for caching
   * inverse operations and therefore need to be used in const methods*
   *.................................................. */
  mutable bool isDecompStale;
  mutable bool areSubmatricesStale;
  mutable Eigen::ColPivHouseholderQR<Eigen::MatrixXd> A_B_decomp;
  mutable Eigen::MatrixXd A_B_cached;
  mutable Eigen::MatrixXd A_N_cached;


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
