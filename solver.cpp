#include <iostream>
#include <cmath>
#include "solver.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Solver(const Point& _borders, const Grid& _grid, double _T, int _TSteps):
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
	double solution = 
		sin( 3 * M_PI * point.x / grid.XSize ) * 
		sin( 3 * M_PI * point.y / grid.YSize ) *
		sin( M_PI * point.z / grid.ZSize) * 
		cos( 4 * (t + M_PI) );
	return solution;
}

void Solver::getAnalyticalSolution(double t, TDArray& result) {
	for (double i = 0; i < grid.XSize(); i++) {
		for (double j = 0; j < grid.YSize(); j++) {
			for (double k = 0; k < grid.ZSize(); k++) {
				Point point = grid.GetPointByIndex(i, j, k);
				result.Value(point) = getAnalyticalSolutionForPoint(point, t);
			}
		}
	}
}

void Solver::PrintOut(double timeLimit, double timeStep) {
	cout << "analytic solution" << endl;

	double time = 0.0;
	while (time < timeLimit) {
		TDArray resultMatrix(grid);
		getAnalyticalSolution(time, resultMatrix);
		cout << "max value is " << resultMatrix.GetMax() << endl;
		cout << "mean value is " << resultMatrix.GetMean() << "\n" << endl;		

		time += timeStep;
	}
}

// Works only for inner points!
double Solver::calcLaplasian(const Point& point, double t, const TDArray& UnValues) {
	double result = 0;
	double middlePart = 2 * UnValues.GetValue({point.x, point.y, point.z});
	result += (UnValues.GetValue({point.x-1, point.y, point.z }) + UnValues.GetValue({point.x+1, point.y, point.z }) - middlePart) / grid.XGridParts;
	result += (UnValues.GetValue({point.x, point.y-1, point.z }) + UnValues.GetValue({point.x, point.y+1, point.z }) - middlePart) / grid.YGridParts;
	result += (UnValues.GetValue({point.x, point.y, point.z-1 }) + UnValues.GetValue({point.x, point.y, point.z+1 }) - middlePart) / grid.ZGridParts;
	return result;
}

double Solver::approximateFunctionInPoint(double laplasian, double UN, double UNMinusOne, double tau) {
	return 2 * Un - UNMinusOne + pow(tau, 2) * laplasian;
}

void Solver::calcUNPlusOne(double time) {

}

//Works correctly if all the data is prepared correctly
void Solver::solve(TDArray& matrix) {

	assert(TSteps > 2);
	//Calculating analytical solution for first two time steps
	getAnalyticalSolution(0, UNMinOne);
	getAnalyticalSolution(tau, UN);

	for(double time = tau * 2; time < T; time += tau) {
		// Step
		calcUNPlusOne(time);

		// rotating
		TDArray* tmp = UNMinOne;
		UNMinOne = UN;
		UN = UNPlusOne;
		UNPlusOne = tmp;
		UNPlusOne.Clean();
	}
}
