#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "AppliedMethods.h"
using namespace std;

//ifstream in("input.txt");
ofstream out("D:\\1\\Documents\\SplineInterpolation2.txt");


int main()
{
	double a = -1, b = 24;
	//double a = 0, b = 2 * 3.14159265358979323846;
	//int n = 10;
	//int n = 25;
	int n = 40;

	double* x = new double[n + 1];
	double* h = new double[n];

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

	SweepMethod(x, coefM_D1T1, n, h, SplineRepresentationType::throughFirstDerivativeType1);
	SweepMethod(x, coefM_D1T2, n, h, SplineRepresentationType::throughFirstDerivativeType2);
	SweepMethod(x, coefM_D1T3, n, h, SplineRepresentationType::throughFirstDerivativeType3);
	SweepMethod(x, coefM_D1T4, n, h, SplineRepresentationType::throughFirstDerivativeType4);
	SweepMethod(x, coefM_D2T1, n, h, SplineRepresentationType::throughSecondDerivativeType1);
	SweepMethod(x, coefM_D2T2, n, h, SplineRepresentationType::throughSecondDerivativeType2);
	SweepMethod(x, coefM_D2T3, n, h, SplineRepresentationType::throughSecondDerivativeType3);
	SweepMethod(x, coefM_D2T4, n, h, SplineRepresentationType::throughSecondDerivativeType4);

	for (double val = a; val <= b; val += 0.001) {

		out << val << ';' << ExactSolution(val) << ';'
			<< BuildingSpline(x, coefM_D1T1, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) << ';'
			<< BuildingSpline(x, coefM_D1T2, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) << ';'
			<< BuildingSpline(x, coefM_D1T3, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) << ';'
			<< BuildingSpline(x, coefM_D1T4, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) << ';'
			<< BuildingSpline(x, coefM_D2T1, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) << ';'
			<< BuildingSpline(x, coefM_D2T2, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) << ';'
			<< BuildingSpline(x, coefM_D2T3, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) << ';'
			<< BuildingSpline(x, coefM_D2T4, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) << ';' << endl;

		approxValue[0] = fabs(BuildingSpline(x, coefM_D1T1, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) - ExactSolution(val));
		approxValue[1] = fabs(BuildingSpline(x, coefM_D1T2, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) - ExactSolution(val));
		approxValue[2] = fabs(BuildingSpline(x, coefM_D1T3, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) - ExactSolution(val));
		approxValue[3] = fabs(BuildingSpline(x, coefM_D1T4, n, h, val, BuildingSplineType::BuildingSplineUsingFirstDerivative) - ExactSolution(val));
		approxValue[4] = fabs(BuildingSpline(x, coefM_D2T1, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) - ExactSolution(val));
		approxValue[5] = fabs(BuildingSpline(x, coefM_D2T2, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) - ExactSolution(val));
		approxValue[6] = fabs(BuildingSpline(x, coefM_D2T3, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) - ExactSolution(val));
		approxValue[7] = fabs(BuildingSpline(x, coefM_D2T4, n, h, val, BuildingSplineType::BuildingSplineUsingSecondDerivative) - ExactSolution(val));

		for (int i = 0; i < 8; i++) {
			if (approxValue[i] > maxApproxValue[i]) {
				maxApproxValue[i] = approxValue[i];
			}
		}
	}
	cout << "Max approx value \n" << endl;
	for (int i = 0; i < 4; i++) {
		cout << "using first deriv, type " << i + 1 << ": " << maxApproxValue[i] << "\t\t using second deriv, type " << i + 1 << ": " << maxApproxValue[i + 4] << endl;
	}

	out.close();

	delete[] x;
	delete[] h;
	delete[] approxValue;
	delete[] maxApproxValue;
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