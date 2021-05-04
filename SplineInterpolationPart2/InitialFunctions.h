#include <cmath>

#pragma region The exact solution and the derivatives of the spline
/// <summary>
/// Точное значение функции f в точке (x,y)
/// </summary>
double ExactSolution(double x, double y) {
	//example 1
	//return 0.1 * (pow(x, 2) - pow(y - 10, 2)) - 2 * x * cos(y) - sin(y);
	
	//example 2
	//return sin(x) * y + x * cos(y) + x + y;

	//example 3
	return -pow(x, 2) - cos(2 * x) * y + sin(2 * x + y);
	
}
/// <summary>
/// Точное значение производной df/dx
///  </summary>
double PartialDerivative10(double x, double y) {
	//example 1
	//return 0.2 * x - 2 * cos(y);
	
	//example 2
	//return y * cos(x) + cos(y) + 1;

	//example 3
	return 2 * sin(2 * x) * y - 2 * x + 2 * cos(2 * x + y);
}
/// <summary>
/// Точное значение производной d^2f/dx/dy
///  </summary>
double MixedPartialDerivative11(double x, double y) {
	//example 1
	//return 2 * sin(y);

	//example 2
	//return cos(x) - sin(y);

	//example 3
	return -2 * sin(2 * x + y) + 2 * sin(2 * x);
}
/// <summary>
/// Точное значение производной df/dy
///  </summary>
double PartialDerivative01(double x, double y) {
	//example 1
	//return 2 * x * sin(y) - 0.2 * y - cos(y) + 2;
	
	//example 2
	//return sin(x) - x * sin(y) + 1;

	//example 3
	return cos(2 * x + y) - cos(2 * x);
}
/// <summary>
/// Точное значение производной d^2f/dx^2
///  </summary>
double PartialDerivative20(double x, double y) {
	//example 1
	//return 0.2;
	
	//example 2
	//return -y * sin(x);

	//example 3
	return 4 * cos(2 * x) * y - 4 * sin(2 * x + y) - 2;
}
/// <summary>
/// Точное значение производной d^4f/dx^2/dy^2
///  </summary>
double MixedPartialDerivative22(double x, double y) {
	//example 1,2
	return 0;

	//example 3
	return 4 * sin(2 * x + y);
}
/// <summary>
/// Точное значение производной d^2f/dy^2
///  </summary>
double PartialDerivative02(double x, double y) {
	//example 1
	//return 2 * x * cos(y) + sin(y) - 0.2;

	//example 2
	//return -x * cos(y);

	//example 3
	return -sin(2 * x + y);
}
#pragma endregion