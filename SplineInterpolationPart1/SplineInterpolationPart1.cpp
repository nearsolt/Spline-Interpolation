#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "AppliedMethods.h"

using namespace std;

//ifstream in("..\\InputInitialConditions.txt");
ofstream out("..\\SplineInterpolation1.txt");

int main(){

#pragma region temp input

	double a = -1, b = 24;
	//double a = 0, b = 2 * 3.14159265358979323846;
	//int n = 10;
	//int n = 25;
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

	double* coefM_D1T1 = new double[n + 1];
	double* coefM_D1T2 = new double[n + 1];
	double* coefM_D1T3 = new double[n + 1];
	double* coefM_D1T4 = new double[n + 1];
	double* coefM_D2T1 = new double[n + 1];
	double* coefM_D2T2 = new double[n + 1];
	double* coefM_D2T3 = new double[n + 1];
	double* coefM_D2T4 = new double[n + 1];

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

	SweepMethod(x, coefM_D1T1, n, h, SplineRepresentationType::throughFirstDerivativeType1);
	SweepMethod(x, coefM_D1T2, n, h, SplineRepresentationType::throughFirstDerivativeType2);
	SweepMethod(x, coefM_D1T3, n, h, SplineRepresentationType::throughFirstDerivativeType3);
	SweepMethod(x, coefM_D1T4, n, h, SplineRepresentationType::throughFirstDerivativeType4);
	SweepMethod(x, coefM_D2T1, n, h, SplineRepresentationType::throughSecondDerivativeType1);
	SweepMethod(x, coefM_D2T2, n, h, SplineRepresentationType::throughSecondDerivativeType2);
	SweepMethod(x, coefM_D2T3, n, h, SplineRepresentationType::throughSecondDerivativeType3);
	SweepMethod(x, coefM_D2T4, n, h, SplineRepresentationType::throughSecondDerivativeType4);

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
		
		spline[0] = BuildingSpline(x, coefM_D1T1, n, h, valOx, BuildingSplineType::buildingSplineUsingFirstDerivative);
		spline[1] = BuildingSpline(x, coefM_D1T2, n, h, valOx, BuildingSplineType::buildingSplineUsingFirstDerivative);
		spline[2] = BuildingSpline(x, coefM_D1T3, n, h, valOx, BuildingSplineType::buildingSplineUsingFirstDerivative);
		spline[3] = BuildingSpline(x, coefM_D1T4, n, h, valOx, BuildingSplineType::buildingSplineUsingFirstDerivative);
		spline[4] = BuildingSpline(x, coefM_D2T1, n, h, valOx, BuildingSplineType::buildingSplineUsingSecondDerivative);
		spline[5] = BuildingSpline(x, coefM_D2T2, n, h, valOx, BuildingSplineType::buildingSplineUsingSecondDerivative);
		spline[6] = BuildingSpline(x, coefM_D2T3, n, h, valOx, BuildingSplineType::buildingSplineUsingSecondDerivative);
		spline[7] = BuildingSpline(x, coefM_D2T4, n, h, valOx, BuildingSplineType::buildingSplineUsingSecondDerivative);

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
	delete[] spline;
	delete[] coefM_D1T1;
	delete[] coefM_D1T2;
	delete[] coefM_D1T3;
	delete[] coefM_D1T4;
	delete[] coefM_D2T1;
	delete[] coefM_D2T2;
	delete[] coefM_D2T3;
	delete[] coefM_D2T4;

	return 0;
}