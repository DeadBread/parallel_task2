#include<cstdlib>
#include "solver.h"

int main(int argc, char *argv[]) {
	Point border = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3])};
	Grid grid = {atoi(argv[4]), atoi(argv[5]), atoi(argv[6])};
	Solver solver(border, grid);
	solver.PrintOut(1, 0.2);
	return 0;
}