#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "AppliedMethods.h"
using namespace std;

//ifstream in("C:\\Users\\1\\source\\repos\\Diplom\\InputInitialConditions3D.txt");
ofstream out("C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt");

/// <summary>
///  Нахождение значений производных в следующей последовательности: m01 -> m10 -> m11
/// </summary>
void FindingDerivatives(double* x, double nX, double* hX, double*y, double nY, double* hY, double** coefM01, double** coefM10, double** coefM11,
			SplineRepresentationType bonCondForOx, SplineRepresentationType bonCondForOy){

	double* tempArrayOx = new double[nX + 1];
	double* tempArrayOy = new double[nY + 1];
	double* tempArrayOxForSecondCycle = new double[nX + 1];
	
	//Найдем m01 (тип краевых условий bonCondForOy) по переменной  y с фиксированной переменной x
	for (int i = 0; i <= nX; i++) {
		for (int j = 0; j <= nY; j++) {
			tempArrayOy[j] = 0;
		}
		SweepMethod(y, tempArrayOy, nY, hY, bonCondForOy, x[i], FixedVariableType::firstCycleFixedVariableX, tempArrayOxForSecondCycle);
		for (int j = 0; j <= nY; j++) {
			coefM01[i][j] = tempArrayOy[j];
		}
	}
	//Найдем m10 (тип краевых условий bonCondForOx) по переменной x с фиксированной переменной y
	for (int j = 0; j <= nY; j++) {
		for (int i = 0; i <= nX; i++) {
			tempArrayOx[i] = 0;
		}
		SweepMethod(x, tempArrayOx, nX, hX, bonCondForOx, y[j], FixedVariableType::firstCycleFixedVariableY, tempArrayOxForSecondCycle);
		for (int i = 0; i <= nX; i++) {
			coefM10[i][j] = tempArrayOx[i];
		}
	}
	//Найдем m11 (тип краевых условий bonCondForOx) по переменной x с фиксированной переменной y
	for (int j = 0; j <= nY; j++) {
		for (int i = 0; i <= nX; i++) {
			tempArrayOx[i] = 0;
			tempArrayOxForSecondCycle[i] = coefM01[i][j];
		}
		SweepMethod(x, tempArrayOx, nX, hX, bonCondForOx, y[j], FixedVariableType::secondCycleFixedVariableY, tempArrayOxForSecondCycle);
		for (int i = 0; i <= nX; i++) {
			coefM11[i][j] = tempArrayOx[i];
		}
	}

	delete[] tempArrayOx;
	delete[] tempArrayOy;
	delete[] tempArrayOxForSecondCycle;
}


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
	double approxValue1 = 0, maxApproxValue1 = 0;

	double** coefM10 = new double* [nX + 1];
	double** coefM01 = new double* [nX + 1];
	double** coefM11 = new double* [nX + 1];

	for (int i = 0; i <= nX; i++) {
		coefM10[i] = new double[nY];
		coefM01[i] = new double[nY];
		coefM11[i] = new double[nY];
	}

	double** coefM101 = new double* [nX + 1];
	double** coefM011 = new double* [nX + 1];
	double** coefM111 = new double* [nX + 1];

	for (int i = 0; i <= nX; i++) {
		coefM101[i] = new double[nY];
		coefM011[i] = new double[nY];
		coefM111[i] = new double[nY];
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



#pragma region Sweep method and finding values m01, m10, m11 
	////Последовательность нахождения значений производных: m01 -> m10 -> m11
	////Найдем m01 (II типа) по переменной  y с фиксированной переменной x
	//for (int i = 0; i <= nX; i++) {
	//	for (int j = 0; j <= nY; j++) {
	//		tempArrayOy[j] = 0;
	//	}
	//	SweepMethod(y, tempArrayOy, nY, hY, SplineRepresentationType::throughFirstDerivativeType2, x[i], FixedVariableType::firstCycleFixedVariableX, tempArrayOxForSecondCycle);
	//	for (int j = 0; j <= nY; j++) {
	//		coefM01[i][j] = tempArrayOy[j];
	//	}
	//}
	////Найдем m10 (II типа) по переменной x с фиксированной переменной y
	//for (int j = 0; j <= nY; j++) {
	//	for (int i = 0; i <= nX; i++) {
	//		tempArrayOx[i] = 0;
	//	}
	//	SweepMethod(x, tempArrayOx, nX, hX, SplineRepresentationType::throughFirstDerivativeType2, y[j], FixedVariableType::firstCycleFixedVariableY, tempArrayOxForSecondCycle);
	//	for (int i = 0; i <= nX; i++) {
	//		coefM10[i][j] = tempArrayOx[i];
	//	}
	//}
	////Найдем m11 (II типа) по переменной x с фиксированной переменной y
	//for (int j = 0; j <= nY; j++) {
	//	for (int i = 0; i <= nX; i++) {
	//		tempArrayOx[i] = 0;
	//		tempArrayOxForSecondCycle[i] = coefM01[i][j];
	//	}
	//	SweepMethod(x, tempArrayOx, nX, hX, SplineRepresentationType::throughFirstDerivativeType2, y[j], FixedVariableType::secondCycleFixedVariableY, tempArrayOxForSecondCycle);
	//	for (int i = 0; i <= nX; i++) {
	//		coefM11[i][j] = tempArrayOx[i];
	//	}
	//}
#pragma endregion
	
	if (!out.is_open()) {
		cout << "Error, invalid output file";
		return 0;
	}

	FindingDerivatives(x, nX, hX, y, nY, hY, coefM01, coefM10, coefM11, SplineRepresentationType::throughFirstDerivativeType2, SplineRepresentationType::throughFirstDerivativeType2);
	FindingDerivatives(x, nX, hX, y, nY, hY, coefM01, coefM10, coefM11, SplineRepresentationType::throughSecondDerivativeType2, SplineRepresentationType::throughSecondDerivativeType2);

	int newNX = (x[nX] - x[0]) * 32; 
	int newNY = (y[nY] - y[0]) * 32;
	
	double newHX = (x[nX] - x[0]) / newNX;
	double newHY = (y[nY] - y[0]) / newNY;
	double valOx, valOy, spline1,spline2;

	for (int i = 0; i <= newNX; i++) {
		valOx = x[0] + newHX * i;
		for (int j = 0; j <= newNY; j++) {
			valOy = y[0] + newHY * j;

			spline1 = BuildingSpline(x, y, coefM10, coefM01, coefM11, nX, nY, hX, hY, valOx, valOy, BuildingSplineType::buildingSplineUsingFirstDerivative);
			spline2 = BuildingSpline(x, y, coefM101, coefM011, coefM111, nX, nY, hX, hY, valOx, valOy, BuildingSplineType::buildingSplineUsingSecondDerivative);

			//out << valOx << ';' << valOy << ';'<< ExactSolution(valOx, valOy) << ';' << spline1 << ';' << spline2 << ';' << endl;
			out << valOx << ';' << valOy << ';'  << spline1 << ';' << ExactSolution(valOx, valOy) << ';' << endl;
			approxValue = fabs(spline1 - ExactSolution(valOx, valOy));
			approxValue1 = fabs(spline2 - ExactSolution(valOx, valOy));
			
			if (approxValue > maxApproxValue) {
				maxApproxValue = approxValue;
			}
			if (approxValue1 > maxApproxValue1) {
				maxApproxValue1 = approxValue1;
			}
		}
	}
	out.close();

	cout<< "Max approx value:" << endl;
	cout << "first deriv: " << maxApproxValue << "\t\t second deriv: " << maxApproxValue1 << endl;


	delete[] x;
	delete[] y;
	delete[] hX;
	delete[] hY;
	//for (int i = 0; i <= nX; i++) {
	//	delete[] coefM10[i];
	//	delete[] coefM01[i];
	//	delete[] coefM11[i];
	//}
	//delete[] coefM10;
	//delete[] coefM01;
	//delete[] coefM11;


	return 0;
}