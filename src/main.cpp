#include <iostream>
#include "lp_solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  LPSolver solver("1 2 3 4 5\n1 2 3 4 5 -6\n9 8 7 6 5 4");

  solver.solve();

  return 0;
}
