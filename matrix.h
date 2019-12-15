#pragma once
#include <stdio.h>
#include <mpi.h>
#include <cstring>

struct Point
{
	double x;
	double y;
	double z;

	// Point operator+(Point& point);
	Point& operator+=(const Point& point) {
		x+=point.x;
		y+=point.y;
		z+=point.z;
		return *this;
	}

	void Print(int rank) const {
		printf("rank = %d; x=%f, y=%f, z=%f\n", rank, x, y, z);
	}
};

class Grid
{
public:
	Grid(int _N, const Point& _shift, const Point& _borders):
		N(_N),
		borders(_borders),
		shift(_shift)
	{}

	int GetN() const { return N; }

	double Lx() const { return borders.x; }
	double Ly() const { return borders.y; }
	double Lz() const { return borders.z; }

	double Xh() const { return borders.x / (N - 1); }
	double Yh() const { return borders.y / (N - 1); }
	double Zh() const { return borders.z / (N - 1); }

	bool IsPointOnBorder(int x, int y, int z) const;

	Point GetPointByIndex(int x, int y, int z) const;
private:
	int N;
	Point borders;
	Point shift;
};

class BorderMatrix {
public:
	explicit BorderMatrix(int x, int y):
		xSize(x), ySize(y) 
	{
		data = new double[GetSize()]();
	}
	~BorderMatrix() {
		delete[] data;
	}

	void Clean() {
		memset(data, 0, GetSize() * sizeof(double));
	}

	int GetSize() const {
		return xSize * ySize; 
	}

	double& Value(int x, int y);

	double GetValue(int x, int y) const {
		return const_cast<BorderMatrix*>(this)->Value(x,y);
	}

	void Send(int to, const MPI_Comm& comm) {
		MPI_Send(data, GetSize(), MPI_DOUBLE, to, 0, comm);
	}

	void Recv(int to, const MPI_Comm& comm) {
		MPI_Status status;
		MPI_Recv(data, GetSize(), MPI_DOUBLE, to, 0, comm, &status);
	}

	void Exchange(int from, int to, const MPI_Comm& comm) {
		MPI_Status status;
		MPI_Sendrecv_replace(data, GetSize(), MPI_DOUBLE, to,
    		0, from, 0, comm, &status);
		if (from == MPI_PROC_NULL) {
			Clean();
			return;
		}
	}

	void Print(int rank);
private:
	int xSize;
	int ySize;

	double* data;
};

class TDArray {
public:
	explicit TDArray(int local_sizes[]);
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

	int XSize() { return xSize; }
	int YSize() { return ySize; }
	int ZSize() { return zSize; }

	double GetAbsMax() const;
	double GetMean() const;

	void Print(int rank) const;

	void Clean();

	// -= operator
	void Subtract(const TDArray& matrix);
private:
	int xSize;
	int ySize;
	int zSize;
	double* data;
};
