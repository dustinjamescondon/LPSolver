#include "lp_solver.hpp"
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <algorithm>

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
  equational_objective_coefficients.resize(num_non_basic_vars + num_basic_vars);
  b_vector.resize(num_basic_vars);

  { // fill the objective coefficient vector
    for(int i = 0; i < first_line_components.size(); ++i) {
      equational_objective_coefficients(i) = std::stod(first_line_components[i]);
    }

    for(int i = 0; i < num_basic_vars; ++i) {
      equational_objective_coefficients(i + num_non_basic_vars) = 0.0;
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
  std::cout << equational_matrix << std::endl;
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

  std::cout << "Basis indices: " <<  basis_indices << std::endl << "non basic indices: " << non_basis_indices << std::endl;
  /*--------------------------------------------------*/
  std::cout << "A_B:\n" << A_B() << std::endl;
  std::cout << "A_N:\n" << A_N() << std::endl;

  std::cout << "Basic objective coefficients:" << x_B() << std::endl;
  std::cout << "Non-basic objective coefficients:" << x_N() << std::endl;
}

// assume the basis indices are sorted?
MatrixXd LPSolver::A_B() const
{
  std::cout << num_basic_vars << std::endl;
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

VectorXd LPSolver::x_N() const
{
  return equational_objective_coefficients(non_basis_indices);
}


VectorXd LPSolver::z_N() const {
  // first solve (A_B)^T v = c_B
  auto v = A_B().transpose().fullPivLu().solve(c_B());

  // then calculate z_N = A_N^T v - c_N
  return A_N().transpose() * v - c_N();
}

double LPSolver::objective_value() {
  // first solve A_B v = b
  auto v = A_B().fullPivLu().solve(b);
  return c_B().dot(v); // TODO make sure this is right
}

VectorXd LPSolver::x_B() const
{
  return equational_objective_coefficients(basis_indices);
}

void LPSolver::pivot(size_t entering, size_t leaving)
{
  for(size_t i = 0; i < basis_indices.size(); ++i) {
    if(basis_indices[i] == leaving) {
      basis_indices[i] = entering;
    }
  }

  for(size_t i = 0; i < non_basis_indices.size(); ++i) {
    if(non_basis_indices[i] == entering){
      non_basis_indices[i] = leaving;
    }
  }

  std::sort(basis_indices.begin(), basis_indices.end(), std::less<unsigned int>());
  std::sort(non_basis_indices.begin(), non_basis_indices.end(), std::less<unsigned int>());
}

void LPSolver::solve()
{
  std::cout << "This is the dictionary:\n" << equational_matrix << std::endl;
  findInitialFeasibleDictionary();
  std::cout << (isUnbounded() ? "unbounded":"bounded");
  std::cout << std::endl;
}

void LPSolver::findInitialFeasibleDictionary()
{
  /* first check if the current dictionary is feasible */
}

/* Returns if the current dictionary is unbounded or not */
bool LPSolver::isUnbounded() const
{
  return false;
}
