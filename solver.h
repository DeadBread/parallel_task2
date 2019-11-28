#pragma once

#include "matrix.h"

class Solver {
public:
	Solver(const Point& _borders, const Grid& _grid);
	void PrintOut(double timeLimit, double timeStep);

private:
	Point borders;
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

	void solve(TDArray& matrix);

	double calcLaplasian(const Point& point, double t, const TDArray& UnValues);
	double approximateFunctionInPoint(double laplasian, double UN, double UNMinusOne, double tau);
	void calcUNPlusOne(double time);
};