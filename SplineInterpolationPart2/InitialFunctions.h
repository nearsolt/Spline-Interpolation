#include <cmath>

#pragma region SplineDerivativesOfTwoVariables
/// <summary>
/// Точное значение функции f в точке (x,y)
/// </summary>
double ExactSolution(double x, double y) {
	return sin(x) * y + x * cos(y) + x + y;
}
/// <summary>
/// Точное значение производной df/dx
///  </summary>
double PartialDerivative10(double x, double y) {
	return y * cos(x) + cos(y) + 1;
}
/// <summary>
/// Точное значение производной d^2f/dx/dy
///  </summary>
double MixedPartialDerivative11(double x, double y) {
	return cos(x) - sin(y);
}
/// <summary>
/// Точное значение производной df/dy
///  </summary>
double PartialDerivative01(double x, double y) {
	return sin(x) - x * sin(y) + 1;
}
/// <summary>
/// Точное значение производной d^2f/dx^2
///  </summary>
double PartialDerivative20(double x, double y) {
	return -y * sin(x);
}
/// <summary>
/// Точное значение производной d^4f/dx^2/dy^2
///  </summary>
double MixedPartialDerivative22(double x, double y) {
	return 0;
}
/// <summary>
/// Точное значение производной d^2f/dy^2
///  </summary>
double PartialDerivative02(double x, double y) {
	return -x * cos(y);
}
#pragma endregion