#include "Calc.h"

Calc::Calc(double dt, double dh) {
	this->dh = dh;
	this->dt = dt;
}

bool Calc::readData(string fileName) {
	ifstream input;
	input.open(fileName, ios::in);
	if (input.is_open())
	{
		input >> N;
		G = new double[N];
		points = new Point[N];
		for (int i = 0; i < N; i++) {
			input >> G[i];
		}
		for (int i = 0; i < N; i++) {
			input >> points[i].x;
			input >> points[i].y;
		}
		return true;
	}
	else return false;
}

void Calc::save(string fileName) {
	if (!output.is_open()) output.open(fileName, ios::out);

	for (int i = 0; i < N - 1; i++) {
		output << points[i].x << ';' << points[i].y << ';';
	}
	output << points[N - 1].x << ';' << points[N - 1].y << endl << endl << endl;
}

void Calc::print() {
	for (int i = 0; i < N - 1; i++) {
		cout << points[i].x << ';' << points[i].y << ';';
	}
	cout << points[N - 1].x << ';' << points[N - 1].y << endl << "Hamiltonian=" << ham() << endl;
}

double Calc::ham() {
	double result = 0;
	try {
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {
				result += G[i] * G[j] * log(pow(points[i].x - points[j].x, 2) + pow(points[i].y - points[j].y, 2));
			}
		}
		result *= -1 / (4 * PI);
	}
	catch (int i) {
		result = ham0;
	}
	return result;
}
double Calc::dHam(int pos, Point* points) {
	double result = 0;
	double h = ham();
	if (pos >= N)
		points[pos - N].y += dh;
	else
		points[pos].x += dh;
	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) {
			result += G[i] * G[j] * log(pow(points[i].x - points[j].x, 2) + pow(points[i].y - points[j].y, 2));
		}
	}
	if (pos >= N)
		points[pos - N].y -= dh;
	else
		points[pos].x -= dh;
	result *= -1 / (4 * PI);
	return (result - h) / dh;
}
//Явный метод Эйлера
void Calc::stepEuler() {
	Point* temp = points;
	Point* newPoints = new Point[N];
	for (int i = 0; i < N; i++) {
		newPoints[i].x = points[i].x + dHam(N + i, points) * dt / G[i];
		newPoints[i].y = points[i].y - dHam(i, points) * dt / G[i];
	}
	points = newPoints;
	t += dt;
	delete(temp);
}
//Метод Рунге-Кутты 4 порядка
void Calc::stepRK() {
	double* k1 = new double[2 * N];
	double* k2 = new double[2 * N];
	double* k3 = new double[2 * N];
	double* k4 = new double[2 * N];
	for (int i = 0; i < N; i++) {
		k1[i] = dt * dHam(N + i, points) / G[i];
		k1[N + i] = -dt * dHam(i, points) / G[i];
	}
	for (int j = 0; j < N; j++) {
		points[j].x += k1[j] / 2;
		points[j].y += k1[N + j] / 2;
	}
	for (int i = 0; i < N; i++) {
		k2[i] = dt * dHam(N + i, points) / G[i];
		k2[N + i] = -dt * dHam(i, points) / G[i];
	}
	for (int j = 0; j < N; j++) {
		points[j].x += -k1[j] / 2 + k2[j] / 2;;
		points[j].y += -k1[N + j] / 2 + k2[N + j] / 2;
	}
	for (int i = 0; i < N; i++) {
		k3[i] = dt * dHam(N + i, points) / G[i];
		k3[N + i] = -dt * dHam(i, points) / G[i];
	}
	for (int j = 0; j < N; j++) {
		points[j].x += -k2[j] / 2 + k3[j];
		points[j].y += -k2[N + j] / 2 + k3[N + j];
	}
	for (int i = 0; i < N; i++) {
		k4[i] = dt * dHam(N + i, points) / G[i];
		k4[N + i] = -dt * dHam(i, points) / G[i];

	}
	for (int j = 0; j < N; j++) {
		points[j].x += -k3[j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
		points[j].y += -k3[N + j] + (k1[N + j] + 2 * k2[N + j] + 2 * k3[N + j] + k4[N + j]) / 6;
	}
	t += dt;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
}

//-----Неявные методы Эйлера (начало)------------------------------------
double Calc::f(int n, Point* points1, Point* points2) {
	double result;
	if (n >= N) {
		result = G[n - N] * (points2[n - N].y - points1[n - N].y) / dt + dHam(n - N, points2);
	}
	else {
		result = G[n] * (points2[n].x - points1[n].x) / dt - dHam(n + N, points2);
	}
	return result;
}
double Calc::df(int n, int k, Point* points1, Point* points2) {
	double result;
	if (k >= N)
		points2[k - N].y += dh;
	else
		points2[k].x += dh;
	result = f(n, points1, points2);
	if (k >= N)
		points2[k - N].y -= dh;
	else
		points2[k].x -= dh;
	return result;
}
//обычный
void Calc::stepImplicitEuler()
{
	Slau m(2 * N);
	Point* points2 = new Point[N];
	for (int i = 0; i < N; i++) {
		points2[i].x = points[i].x;
		points2[i].y = points[i].y;
	}
	do {
		for (int i = 0; i < 2 * N; i++) {
			m.B[i] = f(i, points, points2);
			for (int j = 0; j < 2 * N; j++) {
				m(i, j) = (df(i, j, points, points2) - m.B[i]) / dh;
			}
		}
		m.solve();
		for (int i = 0; i < N; i++) {
			points2[i].x -= m.X[i];
			points2[i].y -= m.X[N + i];
		}
	} while (m.norma() > norm);
	t += dt;
	Point* temp = points;
	points = points2;
	delete temp;
}
//симплектический
void Calc::stepSymplecticEuler()
{
	Slau m(N);
	Point* points2 = new Point[N];
	for (int i = 0; i < N; i++) {
		points2[i].x = points[i].x;
		points2[i].y = points[i].y;
	}
	do {
		for (int i = 0; i < N; i++) {
			m.B[i] = f(i, points, points2);
			for (int j = 0; j < N; j++) {
				m(i, j) = (df(i, j, points, points2) - m.B[i]) / dh;
			}
		}
		m.solve();
		for (int i = 0; i < N; i++) {
			points2[i].x -= m.X[i];
		}
	} while (m.norma() > norm);
	for (int i = 0; i < N; i++) {
		points[i].x = points2[i].x;
		points[i].y -= dHam(i, points2) * dh / G[i];
	}
	t += dt;
	delete[] points2;
}
//-----Неявные методы Эйлера (конец)-------------------------------------
// 
//-----симплектический метод Штермера-Верле (начало)--------------------
double Calc::fV1(int n, Point* points1, Point* points2) {
	return G[n] * (points2[n].x - points1[n].x) / dt - dHam(n + N, points2) / 2;
}
double Calc::dfV1(int n, int k, Point* points1, Point* points2) {
	points2[k].x += dh;
	double result = fV1(n, points1, points2);
	points2[k].x -= dh;
	return result;
}
double Calc::fV2(int n, Point* points1, Point* points2) {
	return G[n] * (points2[n].y - points1[n].y) / dt + dHam(n, points1) / 2 + dHam(n, points2) / 2;
}
double Calc::dfV2(int n, int k, Point* points1, Point* points2) {
	double result;
	points2[k].y += dh;
	result = fV2(n, points1, points2);
	points2[k].y -= dh;
	return result;
}

void Calc::stepSymplecticVerle() {
	Slau m(N);
	Point* points2 = new Point[N];
	for (int i = 0; i < N; i++) {
		points2[i].x = points[i].x;
		points2[i].y = points[i].y;
	}
	do {
		for (int i = 0; i < N; i++) {
			m.B[i] = fV1(i, points, points2);
			for (int j = 0; j < N; j++) {
				m(i, j) = (dfV1(i, j, points, points2) - m.B[i]) / dh;
			}
		}
		m.solve();
		for (int i = 0; i < N; i++) {
			points2[i].x -= m.X[i];
		}
	} while (m.norma() > norm);
	for (int i = 0; i < N; i++) {
		points[i].x = points2[i].x;
	}
	do {
		for (int i = 0; i < N; i++) {
			m.B[i] = fV2(i, points, points2);
			for (int j = 0; j < N; j++) {
				m(i, j) = (dfV2(i, j, points, points2) - m.B[i]) / dh;
			}
		}
		m.solve();
		for (int i = 0; i < N; i++) {
			points2[i].y -= m.X[i];
		}
	} while (m.norma() > norm);
	for (int i = 0; i < N; i++) {
		points[i].y = points2[i].y;
		points[i].x -= dHam(i, points2) * dh / 2;
	}
	t += dt;
	delete[] points2;
}
//-----симплектический метод Штермера-Верле (конец)--------------------