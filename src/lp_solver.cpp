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

  /* Load into the objective function space of the matrix */
  coefficient_matrix = Eigen::MatrixXd(lines.size(), first_line_components.size() + 1);
  for(int i = 0; i < first_line_components.size(); ++i) {
    coefficient_matrix(0, i+1) = std::stod(first_line_components[i]);
  }

  /* Populate the rest of the matrix */
  for(size_t row = 1; row < lines.size(); ++row) {
    vector<string> line_components;
    stringstream streamified_line(lines[row]);
    string component;
    while(getline(streamified_line, component, ' ')) {
      line_components.push_back(component);
    }
    for(int col = 0; col < line_components.size(); ++col) {
      coefficient_matrix(row, col) = std::stod(line_components[col]);
    }
  }
}
