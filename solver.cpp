#include <iostream>
#include <cassert>
#include <cmath>
#include "solver.h"
#include <sstream>

using namespace std;

 #define SIMPLE_BORDERS

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Solver(const Grid& _grid, int padded_local_sizes[], 
	const int* _coords, const int* _dimensions,
	double _T, int _TSteps, const MPI_Comm& _comm, int _rank):
		coords(_coords),
		dimensions(_dimensions),
		grid(_grid), 
		T(_T), 
		TSteps(_TSteps),
		comm(_comm), 
		rank(_rank)
{
	tau = T / (TSteps - 1);
	UNMinOne = new TDArray(padded_local_sizes);
	UN = new TDArray(padded_local_sizes);
	UNPlusOne = new TDArray(padded_local_sizes);
}

double Solver::getAnalyticalSolutionForPoint(const Point& point, double t) {

#ifdef SIMPLE_BORDERS
	//SIMPLE BORDERS
	double t_multiplier = M_PI*sqrt(1/pow(grid.Lx(),2) + 1/pow(grid.Ly(), 2) +1/pow(grid.Lz(), 2));
	double solution = 
		sin( M_PI * point.x / grid.Lx() ) *
		sin( M_PI * point.y / grid.Ly() ) *
		sin( M_PI * point.z / grid.Lz() ) * 
		cos( t_multiplier * t );
#else
	double t_multiplier = M_PI*sqrt(9/pow(grid.Lx(),2) + 4/pow(grid.Ly(), 2) + 4/pow(grid.Lz(), 2));
	double solution = 
		sin( 3 * M_PI * point.x / grid.GetN() ) *
		sin( 2 * M_PI * point.y / grid.GetN() ) *
		sin( 2 * M_PI * point.z / grid.GetN() ) * 
		cos( t_multiplier * t + 4 * M_PI );
#endif
	return solution;
}

void Solver::getAnalyticalSolution(double t, TDArray& result) {
	// cout << "ysize=" << result.YSize() << endl;
#pragma omp parallel for	
	for (int i = 1; i < result.XSize() - 1; i++) {
		for (int j = 1; j < result.XSize() - 1; j++) {
			for (int k = 1; k < result.ZSize() - 1; k++) {
				// cout << "rank " << rank << ", " << i << " " << j << " " << k <<endl;
				Point point = grid.GetPointByIndex(i, j, k);
				// point.Print();
				// printf("%f\n", t);
				result.Value(i, j, k) = getAnalyticalSolutionForPoint(point, t);
				// printf("%f\n", result.Value(i, j, k) );
			}
		}
	}
}

// Works only for inner points!
double Solver::calcLaplasian(int x, int y, int z, double t, const TDArray& UnValues) {
	double result = 0;

	int yRightIndex = y + 1;
	int zRightIndex = z + 1;

	double middlePart = 2*UnValues.GetValue(x, y, z);

	result += (UnValues.GetValue(x-1, y, z ) + UnValues.GetValue(x+1, y, z )- middlePart) / pow(grid.Xh(), 2);
	result += (UnValues.GetValue(x, y-1, z ) + UnValues.GetValue(x, yRightIndex, z ) - middlePart) / pow(grid.Yh(), 2);
	result += (UnValues.GetValue(x, y, z-1 ) + UnValues.GetValue(x, y, zRightIndex ) - middlePart) / pow(grid.Zh(), 2);
	return result;
}

double Solver::approximateFunctionInPoint(double laplasian, double UNValue, double UNMinusOneValue) {
	return 2 * UNValue - UNMinusOneValue + pow(tau, 2) * laplasian;
}

//ineffective, might redo later
void Solver::updateUNBorders() {
	// not passing border elements
	int xs = UNPlusOne->XSize() - 2;
	int ys = UNPlusOne->YSize() - 2;
	int zs = UNPlusOne->ZSize() - 2;

	BorderMatrix XUp(ys, zs);
	BorderMatrix XDown(ys, zs);

	BorderMatrix YUp(xs, zs);
	BorderMatrix YDown(xs, zs);

	BorderMatrix ZUp(xs, ys);
	BorderMatrix ZDown(xs, ys);

	//filling border matrices
#pragma omp parallel for
	for (int j = 1; j < UN->YSize() - 1; ++j) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			XUp.Value(j-1, k-1) = UN->Value(UN->XSize()-2, j, k);
			XDown.Value(j-1, k-1) = UN->Value(1, j, k);
		}
	}

#pragma omp parallel for
	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			YUp.Value(i-1, k-1) = UN->Value(i, UN->YSize()-2, k);
			YDown.Value(i-1, k-1) = UN->Value(i, 1, k);
		}
	}

#pragma omp parallel for
	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int j = 1; j < UN->YSize() - 1; ++j) {
			ZUp.Value(i-1, j-1) = UN->Value(i, j, UN->ZSize()-2);
			ZDown.Value(i-1, j-1) = UN->Value(i, j, 1);
		}
	}

	int from = -1;
	int to = -1;

	stringstream ss;
	ss<<rank << endl;

	MPI_Cart_shift(comm, 0, -1, &from, &to);
	ss << "XDown, from " << from << "to" << to << endl;
    XDown.Exchange(from, to, comm);
    MPI_Cart_shift(comm, 0, 1, &from, &to);
	ss << "XUp, from " << from << " to" << to<< endl;
    XUp.Exchange(from, to, comm);

    MPI_Cart_shift(comm, 1, -1, &from, &to);
	ss << "YDown, from " << from << "to" << to << endl;
    YDown.Exchange(from, to, comm);
    MPI_Cart_shift(comm, 1, 1, &from, &to);
	ss << "YUp, from " << from << " to" << to<< endl;
    YUp.Exchange(from, to, comm);

    MPI_Cart_shift(comm, 2, -1, &from, &to);
	ss << "ZDown, from " << from << "to" << to << endl;
    ZDown.Exchange(from, to, comm);
    MPI_Cart_shift(comm, 2, 1, &from, &to);
	ss << "ZUp, from " << from << " to" << to<< endl;
    ZUp.Exchange(from, to, comm);

	//padding tensor
#pragma omp parallel for
    for (int j = 1; j < UN->YSize() - 1; ++j) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			UN->Value(0, j, k) = XUp.Value(j-1, k-1);
			UN->Value(UN->XSize()-1, j, k) = XDown.Value(j-1, k-1);
		}
	}

#pragma omp parallel for
	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			UN->Value(i, 0, k) = YUp.Value(i-1, k-1);
			UN->Value(i, UN->YSize()-1, k) = YDown.Value(i-1, k-1);
		}
	}

#pragma omp parallel for
	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int j = 1; j < UN->YSize() - 1; ++j) {
			UN->Value(i, j, 0) = ZUp.Value(i-1, j-1);
			UN->Value(i, j, UN->ZSize()-1) = ZDown.Value(i-1, j-1);
		}
	}
}

void Solver::updateBorderConditions(double time) {
	//y borders

#ifndef SIMPLE_BORDERS
	int xs = UNPlusOne->XSize() - 2;
	int ys = UNPlusOne->YSize() - 2;
	int zs = UNPlusOne->ZSize() - 2;

	BorderMatrix YBorder(xs, zs);
	BorderMatrix ZBorder(xs, ys);

	//lower border
	if (coords[1] == 0) {
		//sending second y-slice

#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				YBorder.Value(i-1, k-1) = UN->Value(i, 2, k);
			}
		}
		int dst_coords[3] = {coords[0], dimensions[1] - 1, coords[2]};
		int upper_border_rank = -1;
		MPI_Cart_rank(comm, dst_coords, &upper_border_rank);

		//sending U[:,1,:]
		YBorder.Send(upper_border_rank, comm);

		//recieving I[:,0,:] = U[:,n,:]
		YBorder.Recv(upper_border_rank, comm);


#pragma omp parallel for
		//filling in border
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				UNPlusOne->Value(i, 1, k) = YBorder.Value(i-1, k-1);
			}
		}
	}

	//upper border
	if (coords[1] == dimensions[1] - 1) {
		int dst_coords[3] = {coords[0], 0, coords[2]};
		int lower_border_rank = -1;
		MPI_Cart_rank(comm, dst_coords, &lower_border_rank);

		YBorder.Recv(lower_border_rank, comm);

#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				UN->Value(i, UN->YSize() - 1, k) = YBorder.Value(i-1, k-1);
				// calculation laplasian on the border
				int j = UN->YSize() - 2;
				double laplasian = calcLaplasian(i,j , k, time, *UN);
				UNPlusOne->Value(i,j,k) = approximateFunctionInPoint(laplasian, UN->GetValue(i,j,k), UNMinOne->GetValue(i,j,k) );	
				// updating buffer
				YBorder.Value(i-1, k-1) = UNPlusOne->GetValue(i,j,k);
			}
		}

		YBorder.Send(lower_border_rank, comm);
	}

	//lower border
	if (coords[2] == 0) {

#pragma omp parallel for
		//sending second y-slice
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int j = 1; j < UN->YSize() - 1; ++j) {
				ZBorder.Value(i-1, j-1) = UN->Value(i, j, 2);
			}
		}
		int dst_coords[3] = {coords[0], coords[1], dimensions[2]-1};
		int upper_border_rank = -1;
		MPI_Cart_rank(comm, dst_coords, &upper_border_rank);

		//sending U[:,1,:]
		ZBorder.Send(upper_border_rank, comm);

		//recieving I[:,0,:] = U[:,n,:]
		ZBorder.Recv(upper_border_rank, comm);

#pragma omp parallel for
		//filling in border
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int j = 1; j < UN->YSize() - 1; ++j) {
				UNPlusOne->Value(i, j, 1) = ZBorder.Value(i-1, j-1);
			}
		}
	}

	//upper border
	if (coords[2] == dimensions[2] - 1) {
		int dst_coords[3] = {coords[0],coords[1], 0};
		int lower_border_rank = -1;
		MPI_Cart_rank(comm, dst_coords, &lower_border_rank);

		ZBorder.Recv(lower_border_rank, comm);

#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int j = 1; j < UN->YSize() - 1; ++j) {
				UN->Value(i, j, UN->ZSize() - 1) = ZBorder.Value(i-1, j-1);
				// calculation laplasian on the border
				int k = UN->ZSize() - 2;
				double laplasian = calcLaplasian(i,j,k, time, *UN);
				UNPlusOne->Value(i,j,k) = approximateFunctionInPoint(laplasian, UN->GetValue(i,j,k), UNMinOne->GetValue(i,j,k) );	
				// updating buffer
				ZBorder.Value(i-1, j-1) = UNPlusOne->GetValue(i,j,k);
			}
		}

		ZBorder.Send(lower_border_rank, comm);
	}

#else
	//SIMPLE BORDERS
	if (coords[2] == 0) {
#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int j = 1; j < UN->YSize() - 1; ++j) {
				UNPlusOne->Value(i, j, 1) = 0;
			}
		}
	}

	if (coords[2] == dimensions[2] - 1) {
#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int j = 1; j < UN->YSize() - 1; ++j) {
				UNPlusOne->Value(i, j, UN->ZSize()-2) = 0;
			}
		}
	}

	if (coords[1] == 0) {
#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				UNPlusOne->Value(i, 1, k) = 0;
			}
		}
	}

	if (coords[1] == dimensions[1] - 1 ) {
#pragma omp parallel for
		for (int i = 1; i < UN->XSize() - 1; ++i) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				UNPlusOne->Value(i, UN->YSize()-2, k) = 0;
			}
		}
	}

#endif

	//setting X borders to zero

	if (coords[0] == 0) {
#pragma omp parallel for
		for (int j = 1; j < UN->YSize() - 1; ++j) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				UNPlusOne->Value(1, j, k) = 0;
			}
		}
	}

	if (coords[0] == dimensions[0] - 1 ) {
#pragma omp parallel for
		for (int j = 1; j < UN->YSize() - 1; ++j) {
			for (int k = 1; k < UN->ZSize() - 1; ++k) {
				UNPlusOne->Value(UN->XSize()-2, j, k) = 0;
			}
		}
	}
}


void Solver::calcUNPlusOne(double time) {
#pragma omp parallel for
	for (int i = 1; i < UNPlusOne->XSize() - 1; i++) {
		// For border conditions laplacian is calculated in a specific way
		for (int j = 1; j < UNPlusOne->YSize() - 1; j++) {
			for (int k = 1; k < UNPlusOne->ZSize() - 1; k++) {
				if (grid.IsPointOnBorder(i,j,k)) { 
					UNPlusOne->Value(i,j,k) = 0;
				} else {
					double laplasian = calcLaplasian(i,j,k, time, *UN);
					UNPlusOne->Value(i,j,k) = approximateFunctionInPoint(laplasian, UN->GetValue(i,j,k), UNMinOne->GetValue(i,j,k) );					
				}
			}
		}
	}
}

void Solver::printAndCheck(double time) {
	int sizes[3] = {UNPlusOne->XSize(), UNPlusOne->YSize(), UNPlusOne->ZSize()};
	TDArray anal(sizes);
	getAnalyticalSolution(time, anal);
	anal.Subtract(*UNPlusOne);

	int comm_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    double max_error = anal.GetAbsMax();
    double recv_max_error = -1;
    MPI_Reduce(&max_error, &recv_max_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
    	recv_max_error = max(max_error, recv_max_error);
		cout << "rank" << rank << " time " << time << ", maxError " << max_error << endl;
    }
}

//Works correctly if all the data is prepared correctly
void Solver::Solve() {

    int comm_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	assert(TSteps > 2);
	getAnalyticalSolution(0, *UNMinOne);
	getAnalyticalSolution(tau, *UN);


	for(int i = 2; i < TSteps; i++) {

		double time = tau * i;

		updateUNBorders();

		calcUNPlusOne(time);

#ifndef SIMPLE_BORDERS
		updateBorderConditions(time);
#endif

		// Print and check
		printAndCheck(time);

		// rotating
		TDArray* tmp = UNMinOne;
		UNMinOne = UN;
		UN = UNPlusOne;
		UNPlusOne = tmp;
		UNPlusOne->Clean();
	}
}
