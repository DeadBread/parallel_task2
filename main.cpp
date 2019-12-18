#include<cstdlib>
#include "solver.h"
#include <mpi.h>

void create_shift(int coords[], int dims[], int N, int shift[]) {
	int rank = -1;
   	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (rank == 0)
		printf("dimensions - %d, %d, %d\n", dims[0], dims[1], dims[2]);

	shift[0] = 1.*N/dims[0] * coords[0];
	shift[1]= 1.*N/dims[1] * coords[1];
	shift[2] = 1.*N/dims[2] * coords[2];
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
// 	printf("rank %d comm_size %d\n", rank, comm_size);
	
	int dimensions[3] = {0,0,0};
    int coords[] = {0,0,0};
	MPI_Dims_create(comm_size, 3, dimensions);

	//padded!
    int local_sizes[3] = {N / dimensions[0] + 2, 
    				N/dimensions[1] + 2,
    				N/dimensions[2] + 2};
	// printf("dimensions: %d, %d, %d\n", dimensions[0], dimensions[1], dimensions[2]);
	int periods[3] = {0,0,0};
	MPI_Comm MPI_CART_COMM;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, periods, false, &MPI_CART_COMM);
    MPI_Cart_coords(MPI_CART_COMM, rank, 3, coords);

    int shift[3] = {};
    create_shift(coords, dimensions, N, shift);
    Grid grid(N, shift, border);
    Solver solver(grid, local_sizes, coords, dimensions, 0.002, 11, MPI_CART_COMM, rank);

    solver.Solve();

	TDArray array(local_sizes);

    MPI_Finalize();
	
	
	return 0;
}
