#include <cassert>
#include "matrix.h"

TDArray::TDArray(const Grid& dimensionSizes) :
	xSize(dimensionSizes.XGridParts),
	ySize(dimensionSizes.YGridParts), 
	zSize(dimensionSizes.ZGridParts)
{
	//this shouls fill data with zeros, right?
	data = new double(GetSize());
}

TDArray::~TDArray() {
	delete[] data;
}

double& TDArray::Value(const Point& point) {
	assert(point.x < xSize);
	assert(point.y < ySize);
	assert(point.z < zSize);

	int index = point.x + 
				xSize * point.y + 
				xSize * ySize * point.z;
	return data[index];
}

//TODO: vectorization?
void TDArray::subtract(const TDArray& matrix) {
	for(int i = 0; i < GetSize(); i++) {
		data[i] -= matrix.data[i];
	}
}

double TDArray::GetMax() {
	double max = data[0];
	for(int i = 1; i < GetSize(); i++){
		if (data[i] > max) {
			max = data[i];
		}
	}
	return max;
}

// Attention! Possibility of overflow
double TDArray::GetMean() {
	double sum = 0;
	for(int i = 0; i < GetSize(); i++) {
		sum += data[i];
	}
	return sum / GetSize();	
}
