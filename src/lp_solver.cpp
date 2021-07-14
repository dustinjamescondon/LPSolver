#include "lp_solver.hpp"
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

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
      first_line_components.push_back(component);
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

  std::cout <<"Here's the LP matrix in equational form:\n" << equational_matrix << std::endl;
}

// assume the basis indices are sorted?
// TODO do we need this to be a direct reference of the block in A?
MatrixXd LPSolver::A_B() const
{
  MatrixXd m(num_basic_vars, num_basic_vars);
  for(size_t col = 0; col < num_basic_vars; ++col) {
    m.col(col) = equational_matrix.col(basis_indices(col));
  }
  return m;
}

MatrixXd LPSolver::A_N() const
{
  MatrixXd m(num_basic_vars, num_non_basic_vars);
  for(size_t col = 0; col < num_non_basic_vars; ++col) {
    m.col(col) = equational_matrix.col(non_basis_indices(col));
  }
  return m;
}

// When do we need to calculate this?
VectorXd LPSolver::calcX_B() const
{
  return A_B().fullPivLu().solve(b_vector);
}

VectorXd LPSolver::c_B() const {
  return c_vector(basis_indices);
}

VectorXd LPSolver::c_N() const {
  return c_vector(non_basis_indices);
}

VectorXd LPSolver::calcZ_N() const {
  // first solve (A_B)^T v = c_B
  VectorXd v = A_B().transpose().fullPivLu().solve(c_B());

  // then calculate z_N = A_N^T v - c_N
  return A_N().transpose() * v - c_N();
}

double LPSolver::primalObjectiveValue() const {
  // first solve A_B v = b
  VectorXd v = A_B().fullPivLu().solve(b_vector);
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
}

double LPSolver::dualObjectiveValue() const {
  return c_vector(basis_indices).dot(A_B().fullPivLu().solve(b_vector));
}

// TODO this can't be void: need to return if what we found was optimal or unbounded or something
LPSolver::LPResult LPSolver::dualSolve(Eigen::VectorXd const& obj_coeff_vector) {
  {
    z_vector(basis_indices).fill(0.0);

    double v = A_B().transpose().fullPivLu().solve(c_B());
    z_vector(non_basis_indices) = (A_N().transpose() * v) - c_N();
  }
  if(z_vector(non_basis_indices).minCoeff() < 0.0) {
    // optimal solution to the dual, so return as such
    LPResult r;
    r.isUnbounded = false;
    r.isInfeasible = true;
    return r;
  }

  while(true) {
    size_t leaving_index = chooseLeavingVariable_Dual();
    x_vector(basis_indices) = A_B().fullPivLu().solve(b_vector);
    x_vector(non_basis_indices).fill(0.0);

    if(x_vector(basis_indices).minCoeff() >= 0.0) {
      return dualObjectiveValue();
    }
    // Part 2: choose leaving variable
    size_t leaving_index = chooseLeavingVariable_Dual();

    // Part 3: choose entering variable

    // NOTE NOTE NOTE NOTE Continue here NOTE NOTE NOTE NOTE
  }

}

void LPSolver::solve()
{
  { // Initialize the x vector and check for initial feasibility
    auto x_B = x_vector(basis_indices);
    x_B = calcX_B();
    auto x_N = x_vector(non_basis_indices);
    /*  initialize the x values */
    x_N.fill(0.0);
    // any of the basic variables are negative, then we don't have a feasible dictionary
    if(!isPrimalFeasible()) {
      std::cout << "Initial LP is not primal feasible\n";
      if(isDualFeasible()) {
        // then solve the dual LP
        dualSolve(c_vector);
        /*
        if unbounded
           then primal is infeasible
        else
           optimal solution is optimal solution for primal
        */
      } else {
        // Create
      }

    }
  }

  // If here, we have an initially feasible dictionary, so solve the LP
  while(true) {
    auto z_N = z_vector(non_basis_indices);
    auto z_B = z_vector(basis_indices);
    auto x_B = x_vector(basis_indices);
    auto x_N = x_vector(non_basis_indices);

    // Update the z vector
    z_B.fill(0.0);
    z_N = calcZ_N();

    // If every element of Z is non-negative,
    //   then this is the optimal dictionary
    if(z_N.minCoeff() >= 0.0) {
      std::cout << "optimal\n" << primalObjectiveValue() << std::endl;
      std::cout << x_vector.segment(0, num_non_basic_vars).transpose();
      return;
    }

    // Part 2: choose entering variable (Bland's Rule for now)
    auto entering_index = chooseEnteringVariable();
    std::cout << "Entering variable chosen: " << entering_index << std::endl;

    // Part 3: choose leaving variable
    auto highestIncreaseResult = calcHighestIncrease(entering_index);
    size_t leaving_index = highestIncreaseResult.index;
    double t = highestIncreaseResult.maxIncrease;
    if(highestIncreaseResult.unbounded) {
      std::cout << "unbounded";
    }

    std::cout << "Leaving variable chosen: " << leaving_index << std::endl;

    // Part 4: update for next iteration
    x_B -= t * deltaX(leaving_index);
    x_vector(entering_index) = t;
    pivot(entering_index, leaving_index);
    z_vector(basis_indices).fill(0.0);
    z_vector(non_basis_indices) = calcZ_N();
  }
}


// j is the potential entering variable index
VectorXd LPSolver::deltaX(size_t j)  {
  return A_B().fullPivLu().solve(equational_matrix.col(j));
}

// NOTE Assumes all possible entering vars have been checked for unboundedness, so don't bother
// checking here
// TODO rework how the leaving var is returned
LPSolver::HighestIncreaseResult LPSolver::calcHighestIncrease(unsigned entering_index)  {
  HighestIncreaseResult result;
  VectorXd delta_x = deltaX(entering_index);

  result.unbounded = true; // if we don't find a valid delta_x index, then it's unbounded
  result.maxIncrease = std::numeric_limits<double>::max();
  for(size_t i = 0; i < basis_indices.size(); ++i) {
    size_t basis_index = basis_indices[i];
    if(delta_x[i] > 0.0) {
      double fraction = x_vector[basis_index] / delta_x[i];

      if(fraction < result.maxIncrease) {
        result.maxIncrease = fraction;
        result.index = basis_index;
        result.unbounded = false;
      }
    }
  }

  if(!result.unbounded)  {
    std::cout << "In 'calcHighestIncrease()'; highest increase is: " << result.maxIncrease << std::endl;
  } else {
    std::cout << "In 'calcHighestIncrease()'; unbounded\n";
  }

  return result;
}

LPSolver::HighestIncreaseResult LPSolver::calcHighestIncrease_Dual(unsigned leaving_index) {
  HighestIncreaseResult result;

  VectorXd u_vector(basis_indices.size());
  for(size_t k = 0; k < u_vector.size(); ++k) {
    u_vector(k) = (basis_indices(k) == leaving_index ? 1.0 : 0.0);
  }

  VectorXd delta_z(num_basic_vars + num_non_basic_vars);
  delta_z.fill(0.0);
  VectorXd v = A_B().transpose().fullPivLu().solve(u_vector);

  delta_z(non_basis_indices) = -A_N().transpose() * v;

  result.maxIncrease = std::numeric_limits<double>::max();
  result.unbounded = true;
  for(size_t i = 0; i < non_basis_indices.size(); ++i) {
    size_t non_basis_index = non_basis_indices(i);
    if(delta_z(non_basis_index) > 0.0) {
      result.unbounded = false;
      double fraction = z_vector(non_basis_index) / delta_z(non_basis_index);
      if(fraction < result.maxIncrease) {
         result.maxIncrease = fraction;
         result.index = non_basis_index;
      }
    }
  }

  return result;
}

bool LPSolver::isOptimal()  {
  return calcZ_N().minCoeff() >= 0.0;
}


// NOTE: uses Bland's Rule (lowest index)
size_t LPSolver::chooseEnteringVariable() const {
  VectorXd z_n = z_vector(non_basis_indices);
  for(size_t i = 0; i < non_basis_indices.size(); ++i) {
    if(z_n(i) < 0.0) {
      return non_basis_indices[i];
    }
  }
  return 0; // shouldn't get here
}

// assumes it's not dual optimal, so there will be a leaving var
size_t LPSolver::chooseLeavingVariable_Dual() const {
  VectorXd x_B = x_vector(basis_indices);
  for(size_t i = 0; i < basis_indices.size(); ++i) {
    if(x_B(i) < 0.0) {
      return basis_indices(i);
    }
  }
  return 0;
}



bool LPSolver::isDualFeasible() {
  z_vector(basis_indices).fill(0.0);
  VectorXd v = A_B().transpose().fullPivLu().solve(c_B());
  z_vector(non_basis_indices) = (A_N().transpose() * v) - c_N();

  return (z_vector(non_basis_indices).minCoeff() >= 0.0);
}

bool LPSolver::isPrimalFeasible() const {
  VectorXd x_B = calcX_B();

  // if all the X_B coefficients are non-negative, it's infeasible
  return (x_B.minCoeff() >= 0.0);
}

/* Returns if the basis is unbounded or not
   NOTE assumes the current dictionary isn't optimal, so need to check that before calling this
   method */
bool LPSolver::isUnbounded()
{
  VectorXd z_n = calcZ_N();
  // look through all possible entering variables to see if they can be increased arbitrariy large
  for(size_t i = 0; i < z_n.size(); ++i) {
    // if it's negative, that means the basis variable can be increased, so look at this one
    if(z_n[i] < 0.0) {
      VectorXd deltaX_i = deltaX(non_basis_indices[i]);
      if(deltaX_i.minCoeff()<= 0.0)
        return true;
    }
  }
  return false;
}
