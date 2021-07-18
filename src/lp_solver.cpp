#include "lp_solver.hpp"
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>

//#define DEBUG

const double epsilon = 1.0e-10;

// if a value is smaller than epsilon, then set it to zero
VectorXd zeroify(VectorXd vals) {
  for(double& val: vals) {
    if(std::abs(val) < epsilon)
      val = 0.0;
  }

  return vals;
}

using namespace std;

void LPSolver::updateX() {
  x_vector(basis_indices) = solveA_B_x_equals_b(b_vector);
  x_vector(non_basis_indices).fill(0.0);
}

void LPSolver::updateZ(VectorXd const& obj_coeff_vector) {
    z_vector(basis_indices).fill(0.0);

    VectorXd v = solveA_B_transpose_x_equals_b(obj_coeff_vector(basis_indices));
    z_vector(non_basis_indices) = (A_N().transpose() * v) - obj_coeff_vector(non_basis_indices);
}

LPSolver::LPSolver(const char* filename)
{
  std::ifstream t(filename, std::ifstream::in);
  if (!t.is_open()) {
    std::cout<< "Couldn't open " << filename << std::endl;
    exit(1);
  }
  std::stringstream streamified_text;
  streamified_text << t.rdbuf();

  vector<string> lines;

  /*--------------------------------------------------
   * Parse the lines
   ..................................................*/
  {
    string line;

    while(getline(streamified_text, line, '\n')){
      /*--------------------------------------------------
        replace any tabs with spaces so we only have to
        deal with one type of delimiter
       .................................................. */
      for(char & c : line)
        if(c == '\t') c = ' ';
      /*--------------------------------------------------*/
      lines.push_back(line);
    }
  }
  /*-------------------------------------------------
   * Get rid of any empty lines
   *................................................*/
  for(vector<string>::iterator it = lines.begin(); it != lines.end(); ++it) {
    if(it->empty())
      lines.erase(it);
  }

  /*--------------------------------------------------
   * Count the number of rows and columns
   *..................................................*/
  vector<string> first_line_components;
  {
    stringstream streamified_line(lines[0]);
    string component;
    while(getline(streamified_line, component, ' ')){

      // checking if it's empty will get rid of trailing whitespaces
      if(!component.empty()) {
        first_line_components.push_back(component);
      }
    }
  }

  num_non_basic_vars = first_line_components.size();
  num_basic_vars = lines.size() - 1; // subtracting the objective function line


  /*--------------------------------------------------
   * Allocate the vectors and matrices accordingly
   *..................................................*/
  equational_matrix.resize(num_basic_vars, num_basic_vars + num_non_basic_vars);
  c_vector.resize(num_non_basic_vars + num_basic_vars);
  x_vector.resize(num_non_basic_vars + num_basic_vars);
  z_vector.resize(num_non_basic_vars + num_basic_vars);
  b_vector.resize(num_basic_vars);

  /*--------------------------------------------------
   * Fill the objective coefficient vector
   *..................................................*/
  for(int i = 0; i < first_line_components.size(); ++i) {
    c_vector(i) = std::stod(first_line_components[i]);
  }

  for(int i = 0; i < num_basic_vars; ++i) {
    c_vector(i + num_non_basic_vars) = 0.0;
  }

  /*--------------------------------------------------
   * Populate the equational matrix
   *..................................................*/
  for(size_t row = 1; row < lines.size(); ++row) {
    vector<string> line_components;
    stringstream streamified_line(lines[row]);
    string component;

    while(getline(streamified_line, component, ' ')) {
      if(!component.empty())
        line_components.push_back(component);
    }

    /* Put the constant term into the first column, to get it in dictionary form */
    size_t last_col = num_non_basic_vars;
    b_vector(row-1) = std::stod(line_components[last_col]);

    /* Put the coefficients in the matrix, making sure their signs are inverted
      and in the right place */
    for(size_t col = 0; col < line_components.size() -1; ++col) {
      equational_matrix(row-1, col) = std::stod(line_components[col]);
    }
  }

  // initialize the basis block of the matrix with the identity matrix
  equational_matrix.block(0, num_non_basic_vars, num_basic_vars, num_basic_vars) = Eigen::MatrixXd::Identity(num_basic_vars, num_basic_vars);

  /*--------------------------------------------------
   * Initialize the basis and non-basis lists
   *..................................................*/
  basis_indices.resize(num_basic_vars);
  non_basis_indices.resize(num_non_basic_vars);

  for(int i = 0; i < num_basic_vars; ++i) {
    basis_indices(i) = i + num_non_basic_vars;
  }
  for(int i = 0; i < num_non_basic_vars; ++i) {
    non_basis_indices(i) = i;
  }
  /*--------------------------------------------------*/

  isDecompStale = true;
  areSubmatricesStale = true;
  A_B_cached.resize(num_basic_vars, num_basic_vars);
  A_N_cached.resize(num_basic_vars, num_non_basic_vars);
}

void LPSolver::updateSubmatrices() const {
  for(size_t col = 0; col < num_basic_vars; ++col) {
    A_B_cached.col(col) = equational_matrix.col(basis_indices(col));
  }
  for(size_t col = 0; col < num_non_basic_vars; ++col) {
    A_N_cached.col(col) = equational_matrix.col(non_basis_indices(col));
  }
}

MatrixXd LPSolver::A_B() const
{
  if(areSubmatricesStale) {
    updateSubmatrices();
    areSubmatricesStale = false;
  }

  return A_B_cached;
}

// TODO: i'm repeating my self here
MatrixXd LPSolver::A_N() const
{
  if(areSubmatricesStale) {
    updateSubmatrices();
    areSubmatricesStale = false;
  }
  return A_N_cached;
}

VectorXd LPSolver::c_B() const {
  return c_vector(basis_indices);
}

VectorXd LPSolver::c_N() const {
  return c_vector(non_basis_indices);
}


double LPSolver::primalObjectiveValue() const {
  // first solve A_B v = b
  VectorXd v = solveA_B_x_equals_b(b_vector);
  return c_vector(basis_indices).dot(v);
}

void LPSolver::pivot(size_t entering, size_t leaving)
{
  for(size_t i = 0; i < basis_indices.size(); ++i) {
    if(basis_indices[i] == leaving) {
      basis_indices[i] = entering;
    }
  }

  for(size_t i = 0; i < non_basis_indices.size(); ++i) {
    if(non_basis_indices[i] == entering) {
      non_basis_indices[i] = leaving;
    }
  }

  std::sort(basis_indices.begin(), basis_indices.end(), std::less<unsigned int>());
  std::sort(non_basis_indices.begin(), non_basis_indices.end(), std::less<unsigned int>());

  isDecompStale = true;
  areSubmatricesStale = true;
}

// Assumes there's at least one variable
void LPSolver::printOptimalVariableAssignment() const {
  // if values are very small, just round them to zero for displaying
  auto zeroified_x_vector = zeroify(x_vector);

  // print the first component without a leading space
  printf("%.7g", zeroified_x_vector(0));

  for(double val : zeroified_x_vector.segment(1, num_non_basic_vars - 1)) {
    printf(" %.7g", val);
  }
  printf("\n");
}

double LPSolver::dualObjectiveValue(VectorXd const& obj_coeff_vector) const {
  return obj_coeff_vector(basis_indices).dot(solveA_B_x_equals_b(b_vector));
}

VectorXd LPSolver::solveA_B_transpose_x_equals_b(VectorXd const& b) const {
  if(isDecompStale) {
    A_B_decomp = A_B().colPivHouseholderQr();
    isDecompStale = false;
  }
  return A_B_decomp.transpose().solve(b);
}


VectorXd LPSolver::solveA_B_x_equals_b(VectorXd const& b) const {
  if(isDecompStale) {
    A_B_decomp = A_B().colPivHouseholderQr();
    isDecompStale = false;
  }

  return A_B_decomp.solve(b);
}


LPSolver::LPResult LPSolver::dualSolve(Eigen::VectorXd const& obj_coeff_vector) {
  {
    /*-------------------------------------------------------------
     * Populate z_vector with the current values
     * .............................................................*/
    updateX();
    updateZ(obj_coeff_vector);
    /*-------------------------------------------------------------*/

    // check for initial feasibility
    if(!isDualFeasible()) {
      LPResult r;
      r.state = State::Infeasible;
      return r;
    }
  }

  // If here, it is initially feasible, so solve it
  while(true) {

    /*--------------------------------------------------
     * Part 1: Check for optimality
     *..................................................*/
    if(isDualOptimal()) {
      LPResult r; r.optimal_val = dualObjectiveValue(obj_coeff_vector); r.state = State::Optimal; return r;
    }

    /*--------------------------------------------------
     * Part 2: Choose leaving variable
     *..................................................*/
    size_t leaving_index = chooseDualLeavingVariable_blandsRule();
    VectorXd delta_z = deltaZ(leaving_index);

    /*-------------------------------------------------
     * Part 3: choose entering variable
     *..................................................*/
    auto highestIncreaseResult = calcHighestIncrease_Dual(delta_z);
    if(highestIncreaseResult.unbounded) {
      LPResult r; r.state = State::Unbounded; return r;
    }

    double s = highestIncreaseResult.maxIncrease;
    size_t entering_index = 1234567;//highestIncreaseResult.index;

    for(size_t non_basis_index : non_basis_indices) {
      if(delta_z(non_basis_index) > epsilon) {
        if(std::abs((z_vector(non_basis_index) / delta_z(non_basis_index)) - s) < epsilon) {
          entering_index = non_basis_index;
          break;
        }
      }
    }
    //printf("entering %i    leaving %i\n", entering_index, leaving_index);

    /*--------------------------------------------------
     * Part 4: update for the next iteration
     *.................................................. */
    pivot(entering_index, leaving_index);

    updateX();
    updateZ(obj_coeff_vector);
  }
}

LPSolver::LPResult LPSolver::primalSolve() {
// Initialize the x vector and check for initial feasibility
  updateX();
  updateZ(c_vector);

  if(!isPrimalFeasible()) {
    LPResult r;
    r.state = State::Infeasible;
    return r;
  }

  // If here, we have an initially feasible dictionary, so solve the LP
  while(true) {
    // If every element of Z is non-negative,
    //   then this is the optimal dictionary
    if(isPrimalOptimal()) {
      LPResult r; r.optimal_val = primalObjectiveValue(); r.state = State::Optimal; return r;
    }

    // Part 2: choose entering variable
    auto entering_index = choosePrimalEnteringVariable_blandsRule();

    // Part 3: choose leaving variable
    VectorXd delta_x = deltaX(entering_index);
    auto highestIncreaseResult = calcHighestIncrease(delta_x);

    if(highestIncreaseResult.unbounded) {
      LPResult r;
      r.state = State::Unbounded;
      return r;
    }

    size_t leaving_index = 1234567; //highestIncreaseResult.index;
    double t = highestIncreaseResult.maxIncrease;

    for(size_t basis_index : basis_indices) {
      if(delta_x(basis_index) > epsilon) {
        if(std::abs((x_vector(basis_index)/delta_x(basis_index)) - t) < epsilon) {
          leaving_index = basis_index;
          break;
        }
      }
    }

    //printf("primal: entering: %i   leaving %i\n", entering_index, leaving_index);

    // Part 4: update for next iteration
    pivot(entering_index, leaving_index);

    updateX();
    updateZ(c_vector);
  }
}


void LPSolver::solve()
{
  updateX();
  updateZ(c_vector);

  // any of the basic variables are negative, then we don't have a feasible dictionary
  if(!isPrimalFeasible()) {
    #ifdef DEBUG
    std::cout << "Initial LP is not primal feasible\n";
    #endif
    if(isDualFeasible()) {
      #ifdef DEBUG
      std::cout << "But it is dual feasible!\n";
      #endif
      // then solve the dual LP
      auto result = dualSolve(c_vector);
      if(result.state == State::Unbounded) {
        std::cout << "infeasible";
        return;
      }
      else {
        printf("optimal\n%.7g\n", result.optimal_val);
        printOptimalVariableAssignment();
        return;
      }
    }
    // Otherwise, if it isn't dual-feasible
    else {
      #ifdef DEBUG
      std::cout << "And it isn't intially dual-feasible either...\n";
      #endif
      // Then find a feasible basis
      Eigen::VectorXd zero_vector(c_vector.size());
      zero_vector.fill(0.0);
      auto auxResult = dualSolve(zero_vector);

      if(auxResult.state == State::Optimal) {
        #ifdef DEBUG
        std::cout << "Found the optimal of the aux\n";
        #endif
        auto result = primalSolve();
        if(result.state == State::Unbounded) {
          std::cout <<"unbounded\n";
        }
        else { // otherwise it's optimal
          printf("optimal\n%.7g\n", result.optimal_val);
          printOptimalVariableAssignment();
          return;
        }
      }
      else if(auxResult.state == State::Unbounded) {
#ifdef DEBUG
        std::cout << "Oh dear, the aux is unbounded..." << std::endl;
#endif
        std::cout << "infeasible\n";
      }
      else {// if(auxResult.state == State::Unbounded)
        std::cout << "infeasible\n";
      }
      return;
    }
  }
  auto result = primalSolve();
  if(result.state == State::Unbounded) {
    std::cout << "unbounded";
  } else {
    printf("optimal\n%.7g\n", result.optimal_val);
    printOptimalVariableAssignment();
  }
}


// j is the potential entering variable index
VectorXd LPSolver::deltaX(size_t j) const {
  VectorXd delta_x(num_basic_vars + num_non_basic_vars);
  delta_x(basis_indices) = solveA_B_x_equals_b(equational_matrix.col(j));
  delta_x(non_basis_indices).fill(0.0);
  return delta_x;
}

LPSolver::HighestIncreaseResult LPSolver::calcHighestIncrease(VectorXd const& delta_x) const {
  HighestIncreaseResult result;

  result.unbounded = true; // if we don't find a valid delta_x index, then it's unbounded
  result.maxIncrease = std::numeric_limits<double>::max();
  for(auto basis_index : basis_indices) {
    if(delta_x[basis_index] > epsilon) {
      double fraction = x_vector[basis_index] / delta_x[basis_index];

      if(fraction < result.maxIncrease) {
        result.maxIncrease = fraction;
        result.unbounded = false;
      }
    }
  }

  return result;
}

VectorXd LPSolver::deltaZ(size_t leaving_index) const {

  VectorXd u_vector(basis_indices.size());
  for(size_t k = 0; k < u_vector.size(); ++k) {
    u_vector(k) = (basis_indices(k) == leaving_index ? 1.0 : 0.0);
  }

  // make sure the leaving_index parameter is valid
  assert(u_vector.maxCoeff() == 1.0);

  VectorXd delta_z(num_basic_vars + num_non_basic_vars);
  VectorXd v = solveA_B_transpose_x_equals_b(u_vector);

  delta_z(basis_indices).fill(0.0);
  delta_z(non_basis_indices) = -A_N().transpose() * v;

  return delta_z;
}

LPSolver::HighestIncreaseResult LPSolver::calcHighestIncrease_Dual(VectorXd const& delta_z) const {
  HighestIncreaseResult result;

  result.maxIncrease = std::numeric_limits<double>::max();
  result.unbounded = true;
  for(auto non_basis_index : non_basis_indices) {
    if(delta_z(non_basis_index) > epsilon) {
      result.unbounded = false;
      double fraction = z_vector(non_basis_index) / delta_z(non_basis_index);
      if(fraction < result.maxIncrease) {
         result.maxIncrease = fraction;
      }
    }
  }

  return result;
}

// NOTE: uses Bland's Rule (lowest index index)
size_t LPSolver::choosePrimalEnteringVariable_blandsRule() const {
  for(auto non_basis_index : non_basis_indices) {
    if(z_vector(non_basis_index) < -epsilon) {
      return non_basis_index;
    }
  }
  assert(false);
  return 0; // shouldn't get here
}

// assumes it's not dual optimal, so there will be a leaving var
size_t LPSolver::chooseDualLeavingVariable_blandsRule() const {
  for(auto basis_index : basis_indices) {
    if(x_vector(basis_index) < -epsilon) {
      return basis_index;
    }
  }
  assert(false);
  return 0; // shouldn't get here
}

// assumes it's not dual optimal, so there will be a leaving var
size_t LPSolver::chooseDualLeavingVariable_largestCoeff() const {
  double max = 0.0;
  size_t max_index = 0;
  for(auto basis_index : basis_indices) {
    if(x_vector(basis_index) < 0.0) {
      if(-x_vector(basis_index) > max) {
        max = -x_vector(basis_index);
        max_index = basis_index;
      }
    }
  }
  return max_index;
}

bool LPSolver::isDualFeasible() const {
  return (z_vector(non_basis_indices).minCoeff() > -epsilon);
}

// NOTE: assumes the x_vector is up to date
bool LPSolver::isPrimalFeasible() const {
  // if all the X_B coefficients are non-negative, it's infeasible
  return (x_vector(basis_indices).minCoeff() > -epsilon);
}

bool LPSolver::isPrimalOptimal() const {
  return (z_vector(non_basis_indices).minCoeff() > -epsilon);
}

bool LPSolver::isDualOptimal() const {
  return (x_vector(basis_indices).minCoeff() > -epsilon);
}
