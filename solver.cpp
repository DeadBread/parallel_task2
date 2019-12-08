#include <iostream>
#include <cassert>
#include <cmath>
#include "solver.h"

using namespace std;

#define SIMPLE_BORDERS

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Solver(const Grid& _grid, int padded_local_sizes[], double _T, int _TSteps, const MPI_Comm& _comm, int _rank):
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
	double t_multiplier = M_PI*sqrt(9/pow(grid.Lx(),2) + 9/pow(grid.Ly(), 2) +1/pow(grid.Lz(), 2));
	//SIMPLE BORDERS
	double solution = 
		sin( M_PI * point.x / grid.Lx() ) *
		sin( M_PI * point.y / grid.Ly() ) *
		sin( M_PI * point.z / grid.Lz() ) * 
		cos( t_multiplier * t );
#else
	double solution = 
		sin( 3 * M_PI * point.x / grid.GetN() ) *
		sin( 3 * M_PI * point.y / grid.GetN() ) *
		sin( M_PI * point.z / grid.GetN() ) * 
		cos( 4 * (t + M_PI) );
#endif
	return solution;
}

void Solver::getAnalyticalSolution(double t, TDArray& result) {
	// cout << "ysize=" << result.YSize() << endl;
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

#ifdef SIMPLE_BORDERS
	//SIMPLE BORDERS

	double yRightIndex = y + 1;
	double zRightIndex = z + 1;
#else
	//For border conditions
	double yRightIndex = y % (grid.GetN() - 1) + 1;
	double zRightIndex = z % (grid.GetN() - 1) + 1;
#endif

	// cout << "calculating laplasian for " << x << " " << y << " " << z << endl;

	double middlePart = 2 * UnValues.GetValue(x, y, z);

	result += (UnValues.GetValue(x-1, y, z ) + UnValues.GetValue(x+1, y, z ) - middlePart);// / pow(grid.Xh(), 2);
	// result += (UnValues.GetValue(x, y-1, z ) + UnValues.GetValue(x, yRightIndex, z ) - middlePart) / pow(grid.Yh(), 2);
	// result += (UnValues.GetValue(x, y, z-1 ) + UnValues.GetValue(x, y, zRightIndex ) - middlePart) / pow(grid.Zh(), 2);
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
	for (int j = 1; j < UN->YSize() - 1; ++j) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			XUp.Value(j-1, k-1) = UN->Value(UN->XSize()-2, j, k);
			XDown.Value(j-1, k-1) = UN->Value(1, j, k);
		}
	}

	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			YUp.Value(i-1, k-1) = UN->Value(i, UN->YSize()-2, k);
			YDown.Value(i-1, k-1) = UN->Value(i, 1, k);
		}
	}

	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int j = 1; j < UN->YSize() - 1; ++j) {
			ZUp.Value(i-1, j-1) = UN->Value(i, j, UN->ZSize()-2);
			ZDown.Value(i-1, j-1) = UN->Value(i, j, 1);
		}
	}

	int from = -1;
	int to = -1;

	MPI_Cart_shift(comm, 0, -1, &from, &to);

	// if (from == 0) {
	// 	XDown.Print(from);
	// }

	printf("%d - echangeing\n", rank);
    XDown.Exchange(from, to, comm);
    printf("%d - echanged\n", rank);

 //  	if (from == 0) {
	// 	XDown.Print(to);
	// }

    MPI_Cart_shift(comm, 0, 1, &from, &to);
    XUp.Exchange(from, to, comm);

    MPI_Cart_shift(comm, 1, -1, &from, &to);
    YDown.Exchange(from, to, comm);
    MPI_Cart_shift(comm, 1, 1, &from, &to);
    YUp.Exchange(from, to, comm);

    MPI_Cart_shift(comm, 2, -1, &from, &to);
    ZDown.Exchange(from, to, comm);
    MPI_Cart_shift(comm, 2, 1, &from, &to);
    ZUp.Exchange(from, to, comm);

    // if (rank == )

	//padding tensor
    for (int j = 1; j < UN->YSize() - 1; ++j) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			UN->Value(0, j, k) = XUp.Value(j-1, k-1);
			UN->Value(UN->XSize()-1, j, k) = XDown.Value(j-1, k-1);
		}
	}

	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int k = 1; k < UN->ZSize() - 1; ++k) {
			UN->Value(i, 0, k) = YUp.Value(i-1, k-1);
			UN->Value(i, UN->YSize()-1, k) = YDown.Value(i-1, k-1);
		}
	}

	for (int i = 1; i < UN->XSize() - 1; ++i) {
		for (int j = 1; j < UN->YSize() - 1; ++j) {
			UN->Value(i, j, 0) = ZUp.Value(i-1, j-1);
			UN->Value(i, j, UN->ZSize()-1) = ZDown.Value(i-1, j-1);
		}
	}
}

void Solver::calcUNPlusOne(double time) {
	// UN->Print();

#ifdef SIMPLE_BORDERS

	//SIMPLE BORDERS
	for (int i = 1; i < UNPlusOne->XSize() - 1; i++) {
		// For border conditions laplacian is calculated in a specific way
		for (int j = 1; j < UNPlusOne->YSize() - 1; j++) {
			for (int k = 1; k < UNPlusOne->ZSize() - 1; k++) {
				// cout << "cbefore laplasian " << i << " " << j << " " << k << endl;
				//simple borders

				if (grid.IsPointOnBorder(i,j,k)) { 
					UNPlusOne->Value(i,j,k) = 0;
				} else {
					// Point point = grid.GetPointByIndex(i-1,j-1,k-1);

					double laplasian = calcLaplasian(i,j,k, time, *UN);
					UNPlusOne->Value(i,j,k) = laplasian;//approximateFunctionInPoint(laplasian, UN->GetValue(i,j,k), UNMinOne->GetValue(i,j,k) );					
				}
			}
		}
	}
#else
	// for (int i = 1; i < UNPlusOne->XSize() - 1; i++) {
	// 	// For border conditions laplacian is calculated in a specific way
	// 	for (int j = 1; j < UNPlusOne->YSize(); j++) {
	// 		for (int k = 1; k < UNPlusOne->ZSize(); k++) {
	// 			// cout << "cbefore laplasian " << i << " " << j << " " << k << endl;

	// 			double laplasian = calcLaplasian(i, j, k, time, *UN);
	// 			UNPlusOne->Value(i,j,k) = approximateFunctionInPoint(laplasian, UN->GetValue(i,j,k), UNMinOne->GetValue(i,j,k) );
	// 		}
	// 	}
	// }
	// Border conditions for X will be satisfied already as we set memory to zeros

	//ONLY WRITTER FOR SEQUENTIAL CASE!

	//Border conditions for Y
	// for (double i = 1; i < UNPlusOne->XSize() - 1; i++) {
	// 	for (double k = 1; k < UNPlusOne->ZSize() - 1; k++) {	
	// 		UNPlusOne->Value(i, 0, k) = UNPlusOne->GetValue(i, grid.GetN() - 1, k );
	// 	}
	// }

	// //Border conditions for Z
	// for (double i = 1; i < grid.GetN() - 1; i++) {
	// 	//This time Y axis is complete, so we go all the way through this axis
	// 	for (double j = 1; j < grid.GetN(); j++) {	
	// 		UNPlusOne->Value(i, j, 0) = UNPlusOne->GetValue(i, j, grid.GetN() - 1 );
	// 	}
	// }
#endif
}

void Solver::printAndCheck(double time) {
	int sizes[3] = {UNPlusOne->XSize(), UNPlusOne->YSize(), UNPlusOne->ZSize()};
	TDArray anal(sizes);
	getAnalyticalSolution(time, anal);
	// anal.Subtract(*UNPlusOne);

	int comm_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    double max_error = UNPlusOne->Value(3, 3, 3);
    double recv_max_error = -1;
    MPI_Reduce(&max_error, &recv_max_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
    	recv_max_error = max(max_error, recv_max_error);
		cout << "rank" << rank << " time " << time << ", maxError " << max_error << endl;
    }
}

//Works correctly if all the data is prepared correctly
void Solver::Solve() {

	assert(TSteps > 2);
	//Calculating analytical solution for first two time steps
	getAnalyticalSolution(0, *UNMinOne);
	getAnalyticalSolution(tau, *UN);

	// UNMinOne->Print();
	// return;

	// printf("tau=%f\n", tau);

	// cout << "starting approximation";

	for(int i = 2; i < TSteps; i++) {

		double time = tau * i;
		// Step
		calcUNPlusOne(time);

		updateUNBorders();

		// UNPlusOne->Print(rank);

		// return;

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
