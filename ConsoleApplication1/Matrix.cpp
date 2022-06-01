#include "Matrix.h"
#include <iostream>

Matrix::Matrix(int inSize)
{	
	size = inSize;
	S.resize(size);
	for (int i = 0; i < size; i++) {
		S[i].resize(size,0);			
	}
}
Matrix::~Matrix() {
	//delete[] M;
}
double& Matrix::operator()(int row, int col) const
{	
	return (double&)S[row][col];
}
Matrix Matrix::operator*(const double& f)
{	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			(*this)(i, j) *= f;
		}
	}
	return *this;
}
double Matrix::Determinant() const 
{
	double determ = 0;
	if (size == 1) {
		return (*this)(0, 0);
	}
	for (int i = 0; i < size; i++) {
		double a = (*this)(0, i) * (i % 2 ? 1 : -1);
		determ += a * Minor(0, i).Determinant();
	}
	return determ;
}
Matrix Matrix::Reverse()
{
	Matrix newMatrix(size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			newMatrix(j,i) = Minor(i, j).Determinant() * ((i + j) % 2 ? 1 : -1);
		}
	}
	return newMatrix*(1/Determinant());
}
Matrix Matrix::Minor(int row, int col) const
{
	Matrix newMatrix(size - 1);
	for (int i = 0, in = 0; i < size; i++) {
		if (i != row) {
			for (int j = 0, jn = 0; j < size; j++) {
				if (j != col) {
					newMatrix(in, jn++) = (*this)(i, j);
				}
			}
			in++;
		}
	}
	return newMatrix;
}
void Matrix::Print()
{
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			std::cout << (*this)(i,j)<<' ';
		}
		std::cout << std::endl;
	}
}

Slau::Slau(int inSize):Matrix(inSize)
{
	B.resize(size,0);
	X.resize(size,0);
}
Slau::~Slau() {
	//delete[] M;
	//delete[] B;
}
std::vector<double> Slau::solve() {
	Matrix rev(size);
	rev = this->Reverse();	
	for (int i = 0; i < size; i++) {
		X[i] = 0;
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			X[i] += rev(i,j)*B[j];
		}
	}	
	return X;
}
double Slau::norma() {
	double result = 0;
	for (int i = 0; i < size; i++) {
		result += X[i]>0?X[i]:-X[i];
	}
	return result;
}