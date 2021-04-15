#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "AppliedMethods.h"

using namespace std;

//ifstream in("C:\\Users\\1\\source\\repos\\Diplom\\InputInitialConditions.txt");
ofstream out("C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation1.txt");

int main(){

#pragma region temp input

	double a = -1, b = 24;
	//double a = 0, b = 2 * 3.14159265358979323846;
	int n = 50;
#pragma endregion
	
	//int n;
	//if (!in.is_open()) {
	//	cout << "Error, invalid input file";
	//	return 0;
	//}
	//in >> n;

	double* x = new double[n + 1];
	double* h = new double[n];

	//in >> x[0];
	//for (int i = 1; i <= n; i++) {
	//	in >> x[i];
	//	h[i - 1] = x[i] - x[i - 1];
	//}
	//in.close();
	
	double* approxValue = new double[8]();
	double* maxApproxValue = new double[8]();

	double** coefM = new double* [8];
	for (int i = 0; i < 8; i++) {
		coefM[i] = new double[n + 1];
	}
	double* temp = new double[n + 1];

#pragma region temp initial conditions
	//initial conditions
	x[0] = a;
	x[n] = b;

	double hUn = (b - a) / n;

	for (int i = 0; i < n; i++) {
		h[i] = hUn;
	}

	for (int i = 1; i < n; i++) {
		x[i] = x[i - 1] + h[i - 1];
	}
#pragma endregion

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j <= n; j++) {
			temp[j] = 0;
		}
		SweepMethod(x, temp, n, h, (SplineRepresentationType)i);
		for (int j = 0; j <= n; j++) {
			coefM[i][j] = temp[j];
		}
	}

	if (!out.is_open()) {
		cout << "Error, invalid output file";
		return 0;
	}

	int newN = (x[n] - x[0]) * 1000;
	double newH = (x[n] - x[0]) / newN;
	double valOx;
	double* spline = new double[8];

	for (int i = 0; i <= newN; i++) {
		valOx = x[0] + newH * i;
		for (int j = 0; j < 8; j++) {
			for (int k = 0; k <= n; k++) {
				temp[k] = coefM[j][k];
			}
			if (j < 4) {
				spline[j] = BuildingSpline(x, temp, n, h, valOx, (BuildingSplineType)0);
			}
			else {
				spline[j] = BuildingSpline(x, temp, n, h, valOx, (BuildingSplineType)1);
			}
		}
			
		out << valOx << ';' << ExactSolution(valOx) << ';';
		for (int j = 0; j < 8; j++, out << ';') {
			out << spline[j];
		}
		out << endl;

		for (int j = 0; j < 8; j++) {
			approxValue[j] = fabs(spline[j] - ExactSolution(valOx));
			if (approxValue[j] > maxApproxValue[j]) {
				maxApproxValue[j] = approxValue[j];
			}
		}
	}
	out.close();

	cout << "Max approx value: \n" << endl;
	for (int i = 0; i < 4; i++) {
		cout << "first deriv, type " << i + 1 << ": " << maxApproxValue[i] << "\t\t second deriv, type " << i + 1 << ": " << maxApproxValue[i + 4] << endl;
	}

	delete[] x;
	delete[] h;
	delete[] approxValue;
	delete[] maxApproxValue;

	for (int i = 0; i < 8; i++) {
		delete coefM[i];
	}
	delete[] coefM;
	delete[] temp;
	delete[] spline;

	return 0;
}