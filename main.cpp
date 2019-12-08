#include<cstdlib>
#include "solver.h"
#include <mpi.h>

// Point create_shift(int rank, int dims[], int N) {
// 	Point shift={0,0,0};
// 	double Xh = 1.*N/dims[0];
// 	double Yh = 1.*N/dims[1];
// 	double Zh = 1.*N/dims[2];

// 	if (rank == 0){
// 		printf("dimensions - %d, %d, %d", dims[0], dims[1], dims[2]);
// 	}

// 	while(rank > 0) {
// 		if (rank - dims[1]*dims[2] >= 0) {
// 			rank -= dims[1]*dims[2];
// 			shift.x += Xh;
// 		} else if (rank - dims[2] >= 0) {
// 			rank -= dims[2];
// 			shift.y += Yh;			
// 		} else {
// 			rank -= 1;
// 			shift.z += Zh;
// 		}
// 	}
// 	return shift;
// }


Point create_shift(int coords[], int dims[], int N) {
	Point shift={0,0,0};
	// if (coords[0] == 0)	
	// 	printf("dimensions - %d, %d, %d\n", dims[0], dims[1], dims[2]);

	shift.x = 1.*N/dims[0] * coords[0];
	shift.y = 1.*N/dims[1] * coords[1];
	shift.z = 1.*N/dims[2] * coords[2];
	return shift;
}


int main(int argc, char *argv[]) {
	Point border = {atof(argv[1]), atof(argv[2]), atof(argv[3])};
    int N = atoi(argv[4]);
    // border.Print();

	int comm_size = -1;
	int rank = -1;

	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int dimensions[3] = {0,0,0};
    int coords[] = {0, 0, 0};
	MPI_Dims_create(comm_size, 3, dimensions);

	//padded!
    int local_sizes[3] = {N / dimensions[0] + 2, 
    				N/dimensions[1] + 2,
    				N/dimensions[2] + 2};
	// printf("dimensions: %d, %d, %d\n", dimensions[0], dimensions[1], dimensions[2]);
	int periods[3] = {0,0,0};
	MPI_Comm MPI_CART_COMM;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, periods, true, &MPI_CART_COMM);
    MPI_Cart_coords(MPI_CART_COMM, rank, 3, coords);

    Point shift = create_shift(coords, dimensions, N);
    Grid grid(N, local_sizes, shift, border);
    Solver solver(grid, local_sizes, 0.002, 20, MPI_CART_COMM, rank);

    solver.Solve();

	TDArray array(local_sizes);

    MPI_Finalize();
	// Solver solver(grid, 0.002, 20);
	// solver.Solve();
	return 0;
}