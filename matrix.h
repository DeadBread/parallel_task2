#pragma once

struct Point
{
	double x;
	double y;
	double z;
};

class Grid
{
	int XSize() { return XGridParts; }
	int YSize() { return YGridParts; }
	int ZSize() { return ZGridParts; }

	int Xh() { return borders.x / XGridParts; }
	int Yh() { return borders.y / YGridParts; }
	int Zh() { return borders.z / ZGridParts; }

	Point GetPointByIndex(int x, int y, int z) {
		access(x * Hx() <= borders.x);
		access(y * Hy() <= borders.y);
		access(z * Hz() <= borders.z);

		return {x * Hx(), y * Yh(), z*Zh()};
	}
private:
	int XGridParts;
	int YGridParts;
	int ZGridParts;
	Point borders;
};

class TDArray {
public:
	explicit TDArray(const Grid& dimensionSizes);
	~TDArray();
	// matrix data size in doubles
	int GetSize() {
		return xSize * ySize * zSize; 
	}

	// RW access to values
	double& Value(const Point& point);
	// ReadOnle access to values
	double GetValue(const Point& point) {
		return Value(point);
	}

	double GetMax();
	double GetMean();

	void Clean() {
		memset(data, GetSize() * sizeof(double));
	}

	// -= operator
	void subtract(const TDArray& matrix);
private:
	int xSize;
	int ySize;
	int zSize;
	double* data;
};