#include "lp_solver.hpp"
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

LPSolver::LPSolver(const std::string& text)
{
  vector<string> lines;

  /*--------------------------------------------------
   * Parse the lines
   ..................................................*/
  {
    stringstream streamified_text(text);
    string line;

    while(getline(streamified_text, line, '\n')){
      lines.push_back(line);
    }
  }
  /*--------------------------------------------------*/

  /* Get rid of any empty lines */
  for(vector<string>::iterator it = lines.begin(); it != lines.end(); ++it) {
    if(it->empty())
      lines.erase(it);
  }

  /* Count the number of rows and columns */
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

  equational_matrix.resize(num_basic_vars, num_basic_vars + num_non_basic_vars);
  c_vector.resize(num_non_basic_vars + num_basic_vars);
  x_vector.resize(num_non_basic_vars + num_basic_vars);
  z_vector.resize(num_non_basic_vars + num_basic_vars);
  b_vector.resize(num_basic_vars);

  { // fill the objective coefficient vector
    for(int i = 0; i < first_line_components.size(); ++i) {
      c_vector(i) = std::stod(first_line_components[i]);
    }

    for(int i = 0; i < num_basic_vars; ++i) {
      c_vector(i + num_non_basic_vars) = 0.0;
    }
  }

  /* Populate the equational matrix */
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

  std::cout << "number of basic vars: " << num_basic_vars << std::endl;
  std::cout << "number of non-basic vars: " << num_non_basic_vars << std::endl;
  equational_matrix.block(0, num_non_basic_vars, num_basic_vars, num_basic_vars) = Eigen::MatrixXd::Identity(num_basic_vars, num_basic_vars);
  std::cout << "this is the full equational matrix:\n" << equational_matrix << std::endl;

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

  std::cout << "Basis indices: " <<  basis_indices << std::endl << "non basic indices: " << non_basis_indices.transpose() << std::endl;
  /*--------------------------------------------------*/
  std::cout << "A_B:\n" << A_B() << std::endl;
  std::cout << "A_N:\n" << A_N() << std::endl;

  std::cout << "Basic objective coefficients:" << x_B() << std::endl;
  std::cout << "Non-basic objective coefficients:" << x_N() << std::endl;
}

// assume the basis indices are sorted?
// TODO do we need this to be a direct reference of the block in A?
MatrixXd LPSolver::A_B()
{
  MatrixXd m(num_basic_vars, num_basic_vars);
  for(size_t col = 0; col < num_basic_vars; ++col) {
    m.col(col) = equational_matrix.col(basis_indices(col));
  }
  return m;
}

MatrixXd LPSolver::A_N()
{
  MatrixXd m(num_basic_vars, num_non_basic_vars);
  for(size_t col = 0; col < num_non_basic_vars; ++col) {
    m.col(col) = equational_matrix.col(non_basis_indices(col));
  }
  return m;
}

VectorXd LPSolver::x_N()
{
  return x_vector(non_basis_indices);
}

VectorXd LPSolver::x_B() {
  return x_vector(basis_indices);
}

// When do we need to calculate this?
VectorXd LPSolver::calcX_B()
{
  return A_B().fullPivLu().solve(b_vector);
}

VectorXd LPSolver::c_B() {
  return c_vector(basis_indices);
}

VectorXd LPSolver::c_N() {
  return c_vector(non_basis_indices);
}

VectorXd LPSolver::z_B() {
  return z_vector(basis_indices);
}

VectorXd LPSolver::z_N() {
  return z_vector(non_basis_indices);
}

VectorXd LPSolver::calcZ_N() {
  // first solve (A_B)^T v = c_B
  VectorXd v = A_B().transpose().fullPivLu().solve(c_B());

  // then calculate z_N = A_N^T v - c_N
  return A_N().transpose() * v - c_N();
}

double LPSolver::objective_value() {
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

double LPSolver::solve()
{
  {
    auto x_B = x_vector(basis_indices);
    auto x_N = x_vector(non_basis_indices);
    /*  initialize the x values */
    x_B = calcX_B();
    x_N.fill(0.0);
    // any of the basic variables are negative, then we don't have a feasible dictionary
    if(x_B.minCoeff() < 0.0) {
      return -1.0;
    }
  }

  // if here, we have an initially feasible dictionary, so solve the LP
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
      std::cout << "Hey we found the optimal one!\n" << objective_value() << std::endl;
      std::cout << "The x vector is: " << x_vector.transpose() << std::endl;
      std::cout << "The c vector is: " << c_vector.transpose() << std::endl;
      return objective_value();
    }

    // Part 2: choose entering variable (Bland's Rule for now)
    auto entering_index = chooseEnteringVariable();
    std::cout << "Entering variable chosen: " << entering_index << std::endl;

    // Part 3: choose leaving variable
    unsigned int leaving_index;
    double t = calcHighestIncrease(entering_index, leaving_index);
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
double LPSolver::calcHighestIncrease(unsigned entering_index, unsigned& leaving_index_out)  {
  VectorXd delta_x = deltaX(entering_index);

  double min = std::numeric_limits<double>::max();
  for(size_t i = 0; i < basis_indices.size(); ++i) {
    size_t basis_index = basis_indices[i];
    if(delta_x[i] > 0.0) {
      auto result = x_vector[basis_index] / delta_x[i];
      if(result < min) {
        min = result;
        leaving_index_out = basis_index;
      }
    }
  }
  std::cout << "In 'calcHighestIncrease()'; highest increase is: " << min << std::endl;
  return min; // min is actually the max amount we can increase the given "entering" variable
}

bool LPSolver::isOptimal()  {
  return calcZ_N().minCoeff() >= 0.0;
}

void LPSolver::findInitialFeasibleDictionary()
{
  /* first check if the current dictionary is feasible */
}

// NOTE: uses Bland's Rule (lowest index)
size_t LPSolver::chooseEnteringVariable()  {
  VectorXd z_n = z_N();
  for(size_t i = 0; i < non_basis_indices.size(); ++i) {
    if(z_n(i) < 0.0) {
      return non_basis_indices[i];
    }
  }
  return 0; // shouldn't get here
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
