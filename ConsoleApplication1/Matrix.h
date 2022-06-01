#pragma once
#include <vector>
class Matrix
{
public:
	std::vector < std::vector<double> > S;	
	int size;
	Matrix(int inSize);
	~Matrix();
	double& operator	()	(int row, int col) const;
	Matrix	operator	*	(const double& f);	
	double Determinant() const;
	Matrix Reverse();
	Matrix Minor(int row, int col) const;
	void Print();
};

class Slau:public Matrix {
public:
	std::vector<double> B;
	std::vector<double> X;
	Slau(int inSize);
	~Slau();
	std::vector<double> solve();
	double norma();
};