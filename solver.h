#pragma once

#include "matrix.h"

class Solver {
public:
	Solver(const Grid& _grid, int local_sizes[],
		const int* coords, const int* dimensions,
		double _T, int _TSteps, const MPI_Comm& _comm, int rank);
	
	void Solve();

private:
	Grid grid;

//owns
	TDArray* UNMinOne;
	TDArray* UN;
	TDArray* UNPlusOne;

//external 
	const int* coords;
	const int* dimensions;

	const MPI_Comm& comm;
	int rank;

	double T;
	int TSteps;
	double tau;

	void createGrid();

	void shift(int axis, int direction, int& from, int& to );

	double getAnalyticalSolutionForPoint(const Point& point, double t);
	void getAnalyticalSolution(double t, TDArray& result);

	void printAndCheck(double time);

	void updateUNBorders();
	void updateBorderConditions(double time);

	double calcLaplasian(int x, int y, int z, double t, const TDArray& UnValues);
	double approximateFunctionInPoint(double laplasian, double UN, double UNMinusOne);
	void calcUNPlusOne(double time);
};