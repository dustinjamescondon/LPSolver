##
# LP solver
#
# @file
# @version 0.1

INC=-Iinclude -Imodules/eigen
DEBUG=-g
main: src/main.cpp bin/lp_solver.o
	g++ $(INC) $(DEBUG) src/main.cpp bin/lp_solver.o -o bin/main
bin/lp_solver.o: include/lp_solver.hpp src/lp_solver.cpp
	g++ $(INC) $(DEBUG) -c src/lp_solver.cpp -o bin/lp_solver.o
clean:
	rm ./bin/*


# end
