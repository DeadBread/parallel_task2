#pragma once

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
		XGridParts(x), 
		YGridParts(y),
		ZGridParts(z),
		borders(_borders) 
	{}

	int XSize() const { return XGridParts; }
	int YSize() const { return YGridParts; }
	int ZSize() const { return ZGridParts; }

	double Xh() const  { return borders.x / XGridParts; }
	double Yh() const { return borders.y / YGridParts; }
	double Zh() const { return borders.z / ZGridParts; }

	Point GetPointByIndex(int x, int y, int z) const;
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