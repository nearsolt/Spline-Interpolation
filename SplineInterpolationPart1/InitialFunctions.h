#include <cmath>

#pragma region SplineDerivativesOfSingleVariable
/// <summary>
/// Точное значение функции f в точке x
/// </summary>
double ExactSolution(double x) {
	return cos(3 * x - 7) * log(2 * x + 3) + sin(5 * x - 3);
	//return cos(3*x);
}

/// <summary>
/// Точное значение первой производной
///  </summary>
double FirstDerivative(double x) {
	return 5 * cos(3 - 5 * x) + 3 * log(3 + 2 * x) * sin(7 - 3 * x) + 2 * cos(7 - 3 * x) / (3 + 2 * x);
	//return -3*sin(3*x);
}

/// <summary>
/// Точное значение второй производной
///  </summary>
double SecondDerivative(double x) {
	return 25 * sin(3 - 5 * x) + 12 * sin(7 - 3 * x) / (3 + 2 * x) - cos(7 - 3 * x) * (4 / pow(3 + 2 * x, 2) + 9 * log(3 + 2 * x));
	//return -9*cos(3*x);
}
#pragma endregion