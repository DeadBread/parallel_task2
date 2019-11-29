#pragma once

#include "matrix.h"

class Solver {
public:
	Solver(const Grid& _grid, double _T, int _TSteps);
	
	void Solve();

private:
	Grid grid;

	TDArray* UNMinOne;
	TDArray* UN;
	TDArray* UNPlusOne;

	double T;
	int TSteps;
	double tau;

	void createGrid();

	double getAnalyticalSolutionForPoint(const Point& point, double t);
	void getAnalyticalSolution(double t, TDArray& result);

	void printAndCheck(double time);

	double calcLaplasian(int x, int y, int z, double t, const TDArray& UnValues);
	double approximateFunctionInPoint(double laplasian, double UN, double UNMinusOne);
	void calcUNPlusOne(double time);
};