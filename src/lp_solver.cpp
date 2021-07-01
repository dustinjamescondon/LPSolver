#include "lp_solver.hpp"
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>

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

  size_t num_cols = first_line_components.size() + 1;
  size_t num_rows = lines.size();

  /* Load into the objective function space of the matrix
    NOTE this is putting the first column into dictionary form */
  coefficient_matrix = Eigen::MatrixXd(num_rows, num_cols);

  for(int i = 0; i < first_line_components.size(); ++i) {
    coefficient_matrix(0, i+1) = std::stod(first_line_components[i]);
  }

  /* Populate the rest of the matrix and put into dictionary form at the same time */
  for(size_t row = 1; row < lines.size(); ++row) {
    vector<string> line_components;
    stringstream streamified_line(lines[row]);
    string component;
    while(getline(streamified_line, component, ' ')) {
      line_components.push_back(component);
    }


    /* Put the constant term into the first column, to get it in dictionary form */
    size_t last_col = num_cols - 1;
    coefficient_matrix(row, 0) = std::stod(line_components[last_col]);

    /* Put the coefficients in the matrix, making sure their signs are inverted
      and in the right place */
    for(size_t col = 1; col < line_components.size(); ++col) {
      coefficient_matrix(row, col) = -std::stod(line_components[col-1]);
    }
  }
}


void LPSolver::solve()
{
  std::cout << "This is the dictionary:\n" << coefficient_matrix << std::endl;
  findInitialFeasibleDictionary();
  std::cout << (isUnbounded() ? "unbounded":"bounded");
  std::cout << std::endl;
}

void LPSolver::findInitialFeasibleDictionary()
{
  /* first check if the current dictionary is feasible */
  const auto& constraint_coefficients = coefficient_matrix.block(1,0, coefficient_matrix.rows() - 1, 1);

  if(constraint_coefficients.minCoeff() >= 0.0) {
    std::cout << "initially feasible" << std::endl;
  }
  else {
    /* TODO use some method to find a feasible point */
    std::cout << "not initially feasible" << std::endl;
  }
}


/* Returns if the current dictionary is unbounded or not */
bool LPSolver::isUnbounded() const
{
  /* look at all the potential entering variables */
  const Eigen::MatrixXd& objective_coeffs = coefficient_matrix.block(0,1,1,coefficient_matrix.cols() -1);

  for(size_t i = 0; i < objective_coeffs.cols(); ++i) {
    // if the coeff is positive, it's a potential entering var
    if(objective_coeffs(i) > 0.0) {
      // So check to see if the corresponding column is element-wise non-negative
      const auto& var_coeffs = coefficient_matrix.block(1, i + 1, coefficient_matrix.rows() - 1, 1);
      if(var_coeffs.minCoeff() >= 0.0) {
        //  if it is, then it's unbounded
        return true;
      }
    }
  }

  /* if we're here, then it's bounded */
  return false;
}
