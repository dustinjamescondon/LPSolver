#include <iostream>
#include "lp_solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  LPSolver solver("1 1\n1 0 1\n0 1 1");

  solver.solve();


  return 0;
}
