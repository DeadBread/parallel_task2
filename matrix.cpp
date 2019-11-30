#include <cassert>
#include <cstring>
#include <string>
#include <cstdlib>
#include <stdlib.h>

#include "matrix.h"


Point Grid::GetPointByIndex(int x, int y, int z) const {
	assert(x * Xh() <= borders.x);
	assert(y * Yh() <= borders.y);
	assert(z * Zh() <= borders.z);

	// printf("%d, %f\n", x, Xh());

	Point result = {x * Xh(), y * Yh(), z*Zh()};
	// result.Print();
	return result;
}

const char* get_error_message(int one, int two) {
	char* msg = new char[100];
	sprintf(msg, "%d %d \n", one, two);
	return msg;
}

////////////////////////////////////////////////////////////////////////////////////

TDArray::TDArray(const Grid& dimensionSizes) :
	xSize(dimensionSizes.XSize()),
	ySize(dimensionSizes.YSize()), 
	zSize(dimensionSizes.ZSize())
{
	//this shouls fill data with zeros, right?
	data = new double[GetSize()];
}

void TDArray::Print() const {
	for (int i = 0; i < xSize; i++) {
		// For border conditions laplacian is calculated in a specific way
		for (int j = 0; j < ySize; j++) {
			for (int k = 0; k < zSize; k++) {
				printf("%f ", GetValue(i, j, k));
			}
			printf("\n");
		}
		printf("\n");
	}
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
