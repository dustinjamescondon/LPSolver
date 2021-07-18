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

tests: main
	./bin/main ./tests/simple_initially_feasible.txt
	./bin/main ./tests/simple_initially_infeasible.txt
	./bin/main ./tests/simple_dual_feasible.txt

test_infeasible: main
	./bin/main ./tests/test_LPs/input/vanderbei_exercise2.7.txt

test_hard: main
	./bin/main ./tests/test_LPs/input/netlib_sc50b.txt
test_medium: main
	./bin/main ./tests/test_LPs/input/netlib_afiro.txt
test_loop: main
	./bin/main ./tests/test_LPs/input/netlib_klein1.txt
test_all: main
	./tests/test_vanderbei.sh
	./tests/test_netlib.sh

clean:
	rm ./bin/*


# end
