#include <cassert>
#include <cstring>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>

#include "matrix.h"

using namespace std;

//Local coordinates, blobal point!
Point Grid::GetPointByIndex(int x, int y, int z) const {
	int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// printf("%d, %f\n", x, Xh());
	//transition to the unpadded matrix coords
	// printf("%d before, x=%d, y=%d, z=%d; shift=%f, %f, %f\n", rank, x, y, z, shift.x, shift.y, shift.z);

	x -= 1;
	y -= 1;
	z -= 1;

	// printf("%d after, x=%d, y=%d, z=%d; shift=%f, %f, %f\n", rank, x, y, z, shift.x, shift.y, shift.z);


	//shifting
	x += shift.x;
	y += shift.y;
	z += shift.z;

	// printf("x=%d, y=%d, z=%d; shift=%f, %f, %f\n", x, y, z, shift.x, shift.y, shift.z);

	if (!(x * Xh() <= borders.x)) {
		printf("error! point out of border - x is %d, max=%f\n", x, borders.x);
		assert(false);
	}
	if (!(y * Yh() <= borders.y)) {
		printf("error! point out of border - y is %d, may=%f\n", y, borders.y);
		assert(false);
	}
	if (!(z * Zh() <= borders.z)) {
		printf("error! point out of border - z is %d, maz=%f\n", z, borders.z);
		assert(false);
	}

	// assert(x * Xh() <= borders.x);
	// assert(y * Yh() <= borders.y);
	// assert(z * Zh() <= borders.z);

	Point result = {x * Xh(), y * Yh(), z*Zh()};
	// result += shift;

	// cout << "rank " << rank  << " " << x << " " << y << " " << z << endl;	

	// result.Print(rank);
	return result;
}

bool Grid::IsPointOnBorder(int x, int y, int z) const {
	const Point point = GetPointByIndex(x, y, z);
	return point.x == 0 || point.y == 0 || point.z == 0 ||
			point.x == borders.x ||
			point.y == borders.y ||
			point.z == borders.z;
}
////////////////////////////////////////////////////////////////////////////////////

void BorderMatrix::Print(int rank) {
	stringstream str;
	str << "rank " << rank << endl;
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			str << GetValue(i, j) << " ";
		}
		str << "\n";
	}
	cout << str.str() << endl;
	// printf("%d\n", rank);
}

////////////////////////////////////////////////////////////////////////////////////

TDArray::TDArray(int local_sizes[]) :
	xSize(local_sizes[0]),
	ySize(local_sizes[1]), 
	zSize(local_sizes[2])
{
	//this should fill data with zeros, right?
	data = new double[GetSize()]();
}

void TDArray::Print(int rank) const {
	// char str[3000];
	stringstream str;
	str << "rank " << rank << endl;
	for (int i = 0; i < xSize; i++) {
		// For border conditions laplacian is calculated in a specific way
		for (int j = 0; j < ySize; j++) {
			for (int k = 0; k < zSize; k++) {
				str <<  GetValue(i, j, k) << " ";
			}
			str << "\n";
		}
		str << "\n";
	}
	std::cout << str.str() << endl;;
	// printf("%d\n", rank);
}

TDArray::~TDArray() {
	delete[] data;
}

void TDArray::Clean() {
	memset(data, 0, GetSize() * sizeof(double));
}

double& TDArray::Value(int x, int y, int z) {
	assert(x < xSize);
	assert(y < ySize);
	assert(z < zSize);

	int index = z + 
				zSize * y + 
				ySize * zSize * x;
	return data[index];
}

//TODO: vectorization?
void TDArray::Subtract(const TDArray& matrix) {
	for(int i = 0; i < GetSize(); i++) {
		data[i] -= matrix.data[i];
	}
}

double TDArray::GetAbsMax() const {
	double max = abs(data[0]);
	for(int i = 1; i < GetSize(); i++){
		if (abs(data[i]) > max) {
			max = abs(data[i]);
		}
	}
	return max;
}

// Attention! Possibility of overflow
double TDArray::GetMean() const {
	double sum = 0;
	for(int i = 0; i < GetSize(); i++) {
		sum += data[i];
	}
	return sum / GetSize();	
}
