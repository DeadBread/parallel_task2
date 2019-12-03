#include <cassert>
#include <cstring>
#include <string>
#include <cstdlib>
#include <stdlib.h>

#include "matrix.h"

//Local coordinates, blobal point!
Point Grid::GetPointByIndex(int x, int y, int z) const {
	assert(x * Xh() <= borders.x);
	assert(y * Yh() <= borders.y);
	assert(z * Zh() <= borders.z);

	// printf("%d, %f\n", x, Xh());

	Point result = {x * Xh(), y * Yh(), z*Zh()};
	result += shift;
	// result.Print();
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

TDArray::TDArray(int local_sizes[]) :
	xSize(local_sizes[0]),
	ySize(local_sizes[1]), 
	zSize(local_sizes[2])
{
	//this should fill data with zeros, right?
	data = new double[GetSize()]();
}

void TDArray::Print(int rank) const {
	char str[3000];
	sprintf(str, "rank = %d\n", rank);
	for (int i = 0; i < xSize; i++) {
		// For border conditions laplacian is calculated in a specific way
		for (int j = 0; j < ySize; j++) {
			for (int k = 0; k < zSize; k++) {
				sprintf(str, "%f ", GetValue(i, j, k));
			}
			sprintf(str, "\n");
		}
		sprintf(str, "\n");
	}
	printf("%s\n", str);
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

double TDArray::GetMax() const {
	double max = data[0];
	for(int i = 1; i < GetSize(); i++){
		if (data[i] > max) {
			max = data[i];
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
