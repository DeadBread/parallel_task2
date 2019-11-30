#include<cstdlib>
#include "solver.h"

int main(int argc, char *argv[]) {
	Point border = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3])};
	Grid grid = Grid(atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), border);
	Solver solver(grid, 0.002, 20);
	solver.Solve();
	return 0;
}