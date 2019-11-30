#pragma once
#include <stdio.h>

struct Point
{
	double x;
	double y;
	double z;

	void Print() const {
		printf("x=%f, y=%f, z=%f\n", x, y, z);
	}
};

class Grid
{
public:
	Grid(int x, int y, int z, const Point& _borders):
		XGridPointsCount(x), 
		YGridPointsCount(y),
		ZGridPointsCount(z),
		borders(_borders) 
	{}

	int XSize() const { return XGridPointsCount; }
	int YSize() const { return YGridPointsCount; }
	int ZSize() const { return ZGridPointsCount; }

	double Lx() const { return borders.x; }
	double Ly() const { return borders.y; }
	double Lz() const { return borders.z; }

	double Xh() const { return borders.x / (XGridPointsCount - 1); }
	double Yh() const { return borders.y / (YGridPointsCount - 1); }
	double Zh() const { return borders.z / (ZGridPointsCount - 1); }

	Point GetPointByIndex(int x, int y, int z) const;
private:
	int XGridPointsCount;
	int YGridPointsCount;
	int ZGridPointsCount;
	Point borders;
};

class TDArray {
public:
	explicit TDArray(const Grid& dimensionSizes);
	~TDArray();
	// matrix data size in doubles
	int GetSize() const {
		return xSize * ySize * zSize; 
	}

	// RW access to values
	double& Value(int x, int y, int z);
	// ReadOnle access to values
	double GetValue(int x, int y, int z) const {
		return const_cast<TDArray*>(this)->Value(x,y,z);
	}

	double GetMax() const;
	double GetMean() const;

	void Print() const;

	void Clean();

	// -= operator
	void Subtract(const TDArray& matrix);
private:
	int xSize;
	int ySize;
	int zSize;
	double* data;
};