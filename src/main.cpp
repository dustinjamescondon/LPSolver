#include <iostream>
#include "lp_solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  if(argc == 2) {
    LPSolver solver(argv[1]);
    solver.solve();
  } else {
    std::cout << "Need to provide LP input file using command line argument\n";
  }

  return 0;
}
