#pragma once
#include <iostream>
#include <fstream>
#include <math.h>
#include "Matrix.h"
using namespace std;
struct Point
{
	double x, y;
};

class Calc {
private:
	double dt = 0.1, dh = 0.1, norm = 0.01;
	const double PI = acos(-1.0);
	ofstream output;
	double ham0 = 0;
public:
	int N = 0;
	double t = 0;
	double* G = nullptr;
	Point* points = nullptr;

	Calc(double dt, double dh);
	bool readData(string fileName);
	void save(string fileName);
	void print();
	double ham();
	double dHam(int pos, Point* points);
	//����� ����� ������
	void stepEuler();
	//����� �����-����� 4 �������
	void stepRK();
	//-----������� ������ ������ (������)------------------------------------
	double f(int n, Point* points1, Point* points2);
	double df(int n, int k, Point* points1, Point* points2);
	//�������
	void stepImplicitEuler();
	//���������������
	void stepSymplecticEuler();
	//-----������� ������ ������ (�����)-------------------------------------
	// 
	//-----��������������� ����� ��������-����� (������)--------------------
	double fV1(int n, Point* points1, Point* points2);
	double dfV1(int n, int k, Point* points1, Point* points2);
	double fV2(int n, Point* points1, Point* points2);
	double dfV2(int n, int k, Point* points1, Point* points2);
	void stepSymplecticVerle();
	//-----��������������� ����� ��������-����� (�����)--------------------
};
