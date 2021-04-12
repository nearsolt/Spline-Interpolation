#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "AppliedMethods.h"
using namespace std;

//ifstream in("C:\\Users\\1\\source\\repos\\Diplom\\InputInitialConditions3D.txt");
ofstream out("C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt");



int main() {

#pragma region temp input

	double a = 0, b = 16, c = 0, d = 16;
	int nX = 32, nY = 32;

#pragma endregion

	//int nX, nY;
	//if (!in.is_open()) {
	//	cout << "Error, invalid input file";
	//	return 0;
	//}
	//in >> nX >> nY;

	double* x = new double[nX + 1];
	double* y = new double[nY + 1];

	double* hX = new double[nX];
	double* hY = new double[nY];

	//in >> x[0];
	//for (int i = 1; i <= nX; i++) {
	//	in >> x[i];
	//	hX[i - 1] = x[i] - x[i - 1];
	//}
	//in >> y[0];
	//for (int i = 1; i <= nY; i++) {
	//	in >> y[i];
	//	hY[i - 1] = y[i] - y[i - 1];
	//}
	//in.close();

	double approxValue = 0, maxApproxValue = 0;

	double** coefM10 = new double* [nX + 1];
	double** coefM01 = new double* [nX + 1];
	double** coefM11 = new double* [nX + 1];

	for (int i = 0; i <= nX; i++) {
		coefM10[i] = new double[nY];
		coefM01[i] = new double[nY];
		coefM11[i] = new double[nY];
	}

#pragma region temp initial conditions
	//initial conditions
	x[0] = a;
	x[nX] = b;
	y[0] = c;
	y[nY] = d;

	double hUnX = (b - a) / nX;
	double hUnY = (d - c) / nY;

	for (int i = 0; i < nX; i++) {
		hX[i] = hUnX;
	}
	for (int i = 1; i < nX; i++) {
		x[i] = x[i - 1] + hX[i - 1];
	}

	for (int i = 0; i < nY; i++) {
		hY[i] = hUnY;
	}
	for (int i = 1; i < nY; i++) {
		y[i] = y[i - 1] + hY[i - 1];
	}
#pragma endregion

	double* tempArrayOx = new double[nX + 1];
	double* tempArrayOy = new double[nY + 1];
	double* skip = new double[nX + 1];

#pragma region sweepMetod finding m01,m10,m11
	//найдем m01 (II типа) по переменной  y с фиксированной переменной x
	for (int i = 0; i <= nX; i++) {
		for (int j = 0; j <= nY; j++) {
			tempArrayOy[j] = 0;
		}
		SweepMethod(y, tempArrayOy, nY, hY, SplineRepresentationType::throughFirstDerivativeType2, x[i], FixedVariableType::FixedVariableX, skip);
		for (int j = 0; j <= nY; j++) {
			coefM01[i][j] = tempArrayOy[j];
		}
	}

	//найдем m10 (II типа) по переменной x с фиксированной переменной y
	for (int j = 0; j <= nY; j++) {
		for (int i = 0; i <= nX; i++) {
			tempArrayOx[i] = 0;
		}
		SweepMethod(x, tempArrayOx, nX, hX, SplineRepresentationType::throughFirstDerivativeType2, y[j], FixedVariableType::FixedVariableY, skip);
		for (int i = 0; i <= nX; i++) {
			coefM10[i][j] = tempArrayOx[i];
		}
	}

	//найдем m11 (II типа) по переменной x с фиксированной переменной y
	for (int j = 0; j <= nY; j++) {
		for (int i = 0; i <= nX; i++) {
			tempArrayOx[i] = 0;
			skip[i] = coefM01[i][j];
		}
		SweepMethod(x, tempArrayOx, nX, hX, SplineRepresentationType::throughFirstDerivativeType2, y[j], FixedVariableType::SecondCycle, skip);
		for (int i = 0; i <= nX; i++) {
			coefM11[i][j] = tempArrayOx[i];
		}
	}
#pragma endregion
	

	if (!out.is_open()) {
		cout << "Error, invalid output file";
		return 0;
	}

	int newNX = nX*16, newNY = nY*16;
	
	double newHX = (x[nX] - x[0]) / newNX;
	double newHY = (y[nY] - y[0]) / newNY;
	double valOx, valOy;

	for (int i = 0; i <= newNX; i++) {
		valOx = x[0] + newHX * i;
		for (int j = 0; j <= newNY; j++) {
			valOy = y[0] + newHY * j;
			out << valOx << ';' << valOy << ';' << BuildingSpline(x, y, coefM10, coefM01, coefM11, nX, nY, hX, hY, valOx, valOy, BuildingSplineType::BuildingSplineUsingFirstDerivative)
				<< ';' << ExactSolution(valOx, valOy) << ';' << endl;
			approxValue = fabs(BuildingSpline(x, y, coefM10, coefM01, coefM11, nX, nY, hX, hY, valOx, valOy, BuildingSplineType::BuildingSplineUsingFirstDerivative) - ExactSolution(valOx, valOy));
			if (approxValue > maxApproxValue) {
				maxApproxValue = approxValue;
			}
		}
	}

	out.close();
	cout<< "Max approx value:" << maxApproxValue << endl;

	delete[] x;
	delete[] y;
	delete[] hX;
	delete[] hY;
	/*for (int i = 0; i <= nX; i++) {
		delete[] coefM10[i];
		delete[] coefM11[i];
		delete[] coefM01[i];
	}
	delete[] coefM10;
	delete[] coefM01;
	delete[] coefM11;*/
	delete[] tempArrayOx;
	delete[] tempArrayOy;
	delete[] skip;

	return 0;
}