#include <iostream>
#include <cassert>
#include <cmath>
#include "solver.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Solver(const Grid& _grid, double _T, int _TSteps):
		grid(_grid), 
		T(_T), 
		TSteps(_TSteps)
{
	tau = T / TSteps;
	UNMinOne = new TDArray(grid);
	UN = new TDArray(grid);
	UNPlusOne = new TDArray(grid);
}

double Solver::getAnalyticalSolutionForPoint(const Point& point, double t) {
	point.Print();
	printf("Lx=%d, Ly=%d, Lz=%d\n", grid.XSize(), grid.YSize(), grid.ZSize());
	double solution = 
		sin( 3 * M_PI * point.x / grid.XSize() ) *
		sin( 3 * M_PI * point.y / grid.YSize() ) *
		sin( M_PI * point.z / grid.ZSize() ) * 
		cos( 4 * (t + M_PI) );
	return solution;
	// return point.x + point.y + point.z;
}

void Solver::getAnalyticalSolution(double t, TDArray& result) {
	for (int i = 0; i < grid.XSize(); i++) {
		for (int j = 0; j < grid.YSize(); j++) {
			for (int k = 0; k < grid.ZSize(); k++) {
				Point point = grid.GetPointByIndex(i, j, k);
				// point.Print();
				result.Value(i, j, k) = getAnalyticalSolutionForPoint(point, t);
			}
		}
	}
}

// Works only for inner points!
double Solver::calcLaplasian(int x, int y, int z, double t, const TDArray& UnValues) {
	double result = 0;

	//For border conditions
	double yRightIndex = y % (grid.YSize() - 1) + 1;
	double zRightIndex = z % (grid.ZSize() - 1) + 1;


	// cout << "calculating laplasian for " << x << " " << y << " " << z << endl;

	double middlePart = 2 * UnValues.GetValue(x, y, z);
	result += (UnValues.GetValue(x-1, y, z ) + UnValues.GetValue(x+1, y, z ) - middlePart) / grid.XSize();
	result += (UnValues.GetValue(x, y-1, z ) + UnValues.GetValue(x, yRightIndex, z ) - middlePart) / grid.YSize();
	result += (UnValues.GetValue(x, y, z-1 ) + UnValues.GetValue(x, y, zRightIndex ) - middlePart) / grid.ZSize();
	return result;
}

double Solver::approximateFunctionInPoint(double laplasian, double UNValue, double UNMinusOneValue) {
	return 2 * UNValue - UNMinusOneValue + pow(tau, 2) * laplasian;
}

void Solver::calcUNPlusOne(double time) {
	for (int i = 1; i < grid.XSize() - 1; i++) {
		// For border conditions laplacian is calculated in a specific way
		for (int j = 1; j < grid.YSize(); j++) {
			for (int k = 1; k < grid.ZSize(); k++) {
				// cout << "cbefore laplasian " << i << " " << j << " " << k << endl;


				double laplasian = calcLaplasian(i, j, k, time, *UN);
				UNPlusOne->Value(i,j,k) = approximateFunctionInPoint(laplasian, UN->GetValue(i,j,k), UNMinOne->GetValue(i,j,k) );
			}
		}
	}
	// Border conditions for X will be satisfied already as we set memory to zeros

	//Border conditions for Y
	for (double i = 1; i < grid.XSize() - 1; i++) {
		for (double k = 1; k < grid.ZSize() - 1; k++) {	
			UNPlusOne->Value(i, 0, k) = UNPlusOne->GetValue(i, grid.YSize() - 1, k );
		}
	}

	//Border conditions for Z
	for (double i = 1; i < grid.XSize() - 1; i++) {
		//This time Y axis is complete, so we go all the way through this axis
		for (double j = 1; j < grid.YSize(); j++) {	
			UNPlusOne->Value(i, j, 0) = UNPlusOne->GetValue(i, j, grid.ZSize() - 1 );
		}
	}
}

void Solver::printAndCheck(double time) {
	TDArray anal(grid);
	getAnalyticalSolution(time, anal);
	anal.Subtract(*UNPlusOne);
	cout << "time " << time << ", maxError " << anal.GetMax() << ", mean error " << anal.GetMean() << endl;
}

//Works correctly if all the data is prepared correctly
void Solver::Solve() {

	assert(TSteps > 2);
	//Calculating analytical solution for first two time steps
	getAnalyticalSolution(0, *UNMinOne);
	getAnalyticalSolution(tau, *UN);

	UNMinOne->Print();
	return;
	cout << "starting approximation";

	for(double time = tau * 2; time < T; time += tau) {
		// Step
		calcUNPlusOne(time);

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
