#pragma once

#include "matrix.h"

class Solver {
public:
	Solver(const Grid& _grid, int local_sizes[], double _T, int _TSteps, const MPI_Comm& _comm);
	
	void Solve();

private:
	Grid grid;

	TDArray* UNMinOne;
	TDArray* UN;
	TDArray* UNPlusOne;

	const MPI_Comm& comm;

	double T;
	int TSteps;
	double tau;

	void createGrid();

	double getAnalyticalSolutionForPoint(const Point& point, double t);
	void getAnalyticalSolution(double t, TDArray& result);

	void printAndCheck(double time);

	void updateUNBorders();
	void fillOwnBorders();
	void fillAlianBorders();


	double calcLaplasian(int x, int y, int z, double t, const TDArray& UnValues);
	double approximateFunctionInPoint(double laplasian, double UN, double UNMinusOne);
	void calcUNPlusOne(double time);
};