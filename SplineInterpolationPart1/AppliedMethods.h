#include "InitialFunctions.h"

#pragma region Enumerations
/// <summary>
/// Производная, через которую представляем сплайн, и тип краевых условий
/// </summary>
enum class SplineRepresentationType {
	throughFirstDerivativeType1 = 0,
	throughFirstDerivativeType2 = 1,
	throughFirstDerivativeType3 = 2,
	throughFirstDerivativeType4 = 3,
	throughSecondDerivativeType1 = 4,
	throughSecondDerivativeType2 = 5,
	throughSecondDerivativeType3 = 6,
	throughSecondDerivativeType4 = 7
};

/// <summary>
/// Тип производной, который используется для построение сплайна 
/// </summary>
enum class BuildingSplineType {
	buildingSplineUsingFirstDerivative,
	buildingSplineUsingSecondDerivative
};
#pragma endregion

#pragma region Boundary conditions
/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа I
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType1(double* x, int n, double* h, double** matrix, double* d) {
	matrix[0][0] = 1;												// coef for m_0
	matrix[0][1] = 0;												// coef for m_1
	matrix[n][n - 1] = 0;											// coef for m_N-1
	matrix[n][n] = 1;												// coef for m_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i] / (h[i - 1] + h[i]);				//coef for m_i-1 -  \lambda_i
		matrix[i][i] = 2;											//coef for m_i
		matrix[i][i + 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for m_i+1 -  \mu_i
	}

	d[0] = FirstDerivative(x[0]);
	d[n] = FirstDerivative(x[n]);

	for (int i = 1; i < n; i++) {
		d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1]) - ExactSolution(x[i])) / (h[i - 1] + h[i]) / h[i] + h[i] * (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа II
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType2(double* x, int n, double* h, double** matrix, double* d) {
	matrix[0][0] = 2;												// coef for m_0
	matrix[0][1] = 1;												// coef for m_1
	matrix[n][n - 1] = 1;											// coef for m_N-1
	matrix[n][n] = 2;												// coef for m_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i] / (h[i - 1] + h[i]);				//coef for m_i-1 -  \lambda_i
		matrix[i][i] = 2;											//coef for m_i
		matrix[i][i + 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for m_i+1 -  \mu_i
	}

	d[0] = 3 * (ExactSolution(x[1]) - ExactSolution(x[0])) / h[0] - 0.5 * h[0] * SecondDerivative(x[0]);
	d[n] = 3 * (ExactSolution(x[n]) - ExactSolution(x[n - 1])) / h[n - 1] - 0.5 * h[n - 1] * SecondDerivative(x[n]);

	for (int i = 1; i < n; i++) {
		d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1]) - ExactSolution(x[i])) / (h[i - 1] + h[i]) / h[i] + h[i] * (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа III
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType3(double* x, int n, double* h, double** matrix, double* d) {
	matrix[1][1] = 2;												// coef for m_1
	matrix[1][2] = h[0] / (h[0] + h[1]);							// coef for m_2
	matrix[1][n] = h[1] / (h[0] + h[1]);							// coef for m_N
	matrix[n][1] = h[n - 1] / (h[n - 1] + h[n]);					// coef for m_1
	matrix[n][n - 1] = h[n] / (h[n - 1] + h[n]);					// coef for m_N-1
	matrix[n][n] = 2;												// coef for m_N

	for (int i = 2; i < n; i++) {
		matrix[i][i - 1] = h[i] / (h[i - 1] + h[i]);				//coef for m_i-1 -  \lambda_i
		matrix[i][i] = 2;											//coef for m_i
		matrix[i][i + 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for m_i+1 -  \mu_i
	}

	d[n] = 3 * (h[n - 1] * (ExactSolution(x[1]) - ExactSolution(x[0])) / (h[n - 1] + h[0]) / h[0] + h[0] * (ExactSolution(x[0]) - ExactSolution(x[n - 1])) / (h[n - 1] + h[0]) / h[n - 1]);

	for (int i = 1; i < n; i++) {
		d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1]) - ExactSolution(x[i])) / (h[i - 1] + h[i]) / h[i] + h[i] * (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа IV
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType4(double* x, int n, double* h, double** matrix, double* d) {
	matrix[0][0] = 1;												// coef for m_0
	matrix[0][1] = 1 - pow(h[0] / h[1], 2);							// coef for m_1
	matrix[0][2] = -pow(h[0] / h[1], 2);							// coef for m_2
	matrix[1][1] = 1 + h[0] / h[1];									// coef for m_1
	matrix[1][2] = h[0] / h[1];										// coef for m_2
	matrix[n - 1][n - 2] = h[n - 1] / h[n - 2];						// coef for m_N-2
	matrix[n - 1][n - 1] = 1 + h[n - 1] / h[n - 2];					// coef for m_N-1
	matrix[n][n - 2] = -pow(h[n - 1] / h[n - 2], 2);				// coef for m_N-2
	matrix[n][n - 1] = 1 - pow(h[n - 1] / h[n - 2], 2);				// coef for m_N-1
	matrix[n][n] = 1;												// coef for m_N

	for (int i = 2; i < n - 1; i++) {
		matrix[i][i - 1] = h[i] / (h[i - 1] + h[i]);				//coef for m_i-1 -  \lambda_i
		matrix[i][i] = 2;											//coef for m_i
		matrix[i][i + 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for m_i+1 -  \mu_i
	}

	d[0] = 2 * ((ExactSolution(x[1]) - ExactSolution(x[0])) / h[0] - pow(h[0] / h[1], 2) * (ExactSolution(x[2]) - ExactSolution(x[1])) / h[1]);
	d[1] = (h[0] / (h[0] + h[1]) + 2 * h[0] / h[1]) * (ExactSolution(x[2]) - ExactSolution(x[1])) / h[1] + h[1] * (ExactSolution(x[1]) - ExactSolution(x[0])) / (h[0] + h[1]) / h[0];
	d[n - 1] = h[n - 2] * (ExactSolution(x[n]) - ExactSolution(x[n - 1])) / (h[n - 2] + h[n - 1]) / h[n - 1] + (h[n - 1] / (h[n - 2] + h[n - 1]) + 2 * h[n - 1] / h[n - 2]) * (ExactSolution(x[n - 1]) - ExactSolution(x[n - 2])) / h[n - 2];
	d[n] = 2 * ((ExactSolution(x[n]) - ExactSolution(x[n - 1])) / h[n - 1] - pow(h[n - 1] / h[n - 2], 2) * (ExactSolution(x[n - 1]) - ExactSolution(x[n - 2])) / h[n - 2]);

	for (int i = 2; i < n - 1; i++) {
		d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1]) - ExactSolution(x[i])) / (h[i - 1] + h[i]) / h[i] + h[i] * (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа I
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType1(double* x, int n, double* h, double** matrix, double* d) {
	matrix[0][0] = 2;												// coef for M_0
	matrix[0][1] = 1;												// coef for M_1
	matrix[n][n - 1] = 1;											// coef for M_N-1
	matrix[n][n] = 2;												// coef for M_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for M_i-1 - \mu_i
		matrix[i][i] = 2;											//coef for M_i
		matrix[i][i + 1] = h[i] / (h[i - 1] + h[i]);				//coef for M_i+1 -  \lambda_i
	}

	d[0] = 6 * ((ExactSolution(x[1]) - ExactSolution(x[0])) / h[0] - FirstDerivative(x[0])) / h[0];
	d[n] = 6 * (FirstDerivative(x[n]) - (ExactSolution(x[n]) - ExactSolution(x[n - 1])) / h[n - 1]) / h[n - 1];

	for (int i = 1; i < n; i++) {
		d[i] = 6 * ((ExactSolution(x[i + 1]) - ExactSolution(x[i])) / h[i] - (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа II
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType2(double* x, int n, double* h, double** matrix, double* d) {
	matrix[0][0] = 1;												// coef for M_0
	matrix[0][1] = 0;												// coef for M_1
	matrix[n][n - 1] = 0;											// coef for M_N-1
	matrix[n][n] = 1;												// coef for M_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for M_i-1 - \mu_i
		matrix[i][i] = 2;											//coef for M_i
		matrix[i][i + 1] = h[i] / (h[i - 1] + h[i]);				//coef for M_i+1 -  \lambda_i
	}

	d[0] = SecondDerivative(x[0]);
	d[n] = SecondDerivative(x[n]);

	for (int i = 1; i < n; i++) {
		d[i] = 6 * ((ExactSolution(x[i + 1]) - ExactSolution(x[i])) / h[i] - (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа III
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType3(double* x, int n, double* h, double** matrix, double* d) {
	matrix[1][1] = 2;															// coef for M_1
	matrix[1][2] = h[1] / (h[0] + h[1]);										// coef for M_2
	matrix[1][n] = h[0] / (h[0] + h[1]);										// coef for M_N
	matrix[n][1] = h[n] / (h[n - 1] + h[n]);									// coef for M_1
	matrix[n][n - 1] = h[n - 1] / (h[n - 1] + h[n]);							// coef for M_N-1
	matrix[n][n] = 2;															// coef for M_N

	for (int i = 2; i < n; i++) {
		matrix[i][i - 1] = h[i - 1] / (h[i - 1] + h[i]);						//coef for M_i-1 - \mu_i
		matrix[i][i] = 2;														//coef for M_i
		matrix[i][i + 1] = h[i] / (h[i - 1] + h[i]);							//coef for M_i+1 -  \lambda_i
	}

	d[n] = 6 * ((ExactSolution(x[1]) - ExactSolution(x[0])) / h[0] - (ExactSolution(x[0]) - ExactSolution(x[n - 1])) / h[n - 1]) / (h[n - 1] + h[0]);

	for (int i = 1; i < n; i++) {
		d[i] = 6 * ((ExactSolution(x[i + 1]) - ExactSolution(x[i])) / h[i] - (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа IV
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType4(double* x, int n, double* h, double** matrix, double* d) {
	matrix[0][0] = 1;															// coef for M_0
	matrix[0][1] = -(h[0] + h[1]) / h[1];										// coef for M_1
	matrix[0][2] = h[0] / h[1];													// coef for M_2
	matrix[1][1] = 1 + h[1] / (h[0] + h[1]);									// coef for M_1
	matrix[1][2] = (h[1] - h[0]) / (h[0] + h[1]);								// coef for M_2
	matrix[n - 1][n - 2] = (h[n - 2] - h[n - 1]) / (h[n - 2] + h[n - 1]);		// coef for M_N-2
	matrix[n - 1][n - 1] = 1 + h[n - 2] / (h[n - 2] + h[n - 1]);				// coef for M_N-1
	matrix[n][n - 2] = h[n - 1] / h[n - 2];										// coef for M_N-2
	matrix[n][n - 1] = -(h[n - 2] + h[n - 1]) / h[n - 2];						// coef for M_N-1
	matrix[n][n] = 1;															// coef for M_N

	for (int i = 2; i < n - 1; i++) {
		matrix[i][i - 1] = h[i - 1] / (h[i - 1] + h[i]);						//coef for M_i-1 - \mu_i
		matrix[i][i] = 2;														//coef for M_i
		matrix[i][i + 1] = h[i] / (h[i - 1] + h[i]);							//coef for M_i+1 -  \lambda_i
	}

	d[1] = 6 * h[1] * ((ExactSolution(x[2]) - ExactSolution(x[1])) / h[1] - (ExactSolution(x[1]) - ExactSolution(x[0])) / h[0]) / pow(h[0] + h[1], 2);
	d[n - 1] = 6 * h[n - 2] * ((ExactSolution(x[n]) - ExactSolution(x[n - 1])) / h[n - 1] - (ExactSolution(x[n - 1]) - ExactSolution(x[n - 2])) / h[n - 2]) / pow(h[n - 2] + h[n - 1], 2);

	for (int i = 2; i < n - 1; i++) {
		d[i] = 6 * ((ExactSolution(x[i + 1]) - ExactSolution(x[i])) / h[i] - (ExactSolution(x[i]) - ExactSolution(x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
	}
}
#pragma endregion

#pragma region Sweep method
/// <summary>
/// Метод прогонки
/// </summary>
void SweepMethod(double* x, double* coefM, int n, double* h, SplineRepresentationType bonCond) {
	double* alpha = new double[n]();
	double* beta = new double[n]();

	double* d = new double[n + 1]();
	double** matrix = new double* [n + 1];
	for (int i = 0; i <= n; i++) {
		matrix[i] = new double[n + 1]();
	}

	coefM[0] = ExactSolution(x[0]);
	coefM[n] = ExactSolution(x[n]);

	switch (bonCond) {
		case SplineRepresentationType::throughFirstDerivativeType1:
			BuildingTridiagonalMatrixFirstDerivType1(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughFirstDerivativeType2:
			BuildingTridiagonalMatrixFirstDerivType2(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughFirstDerivativeType3:
			BuildingTridiagonalMatrixFirstDerivType3(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughFirstDerivativeType4:
			BuildingTridiagonalMatrixFirstDerivType4(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughSecondDerivativeType1:
			BuildingTridiagonalMatrixSecondDerivType1(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughSecondDerivativeType2:
			BuildingTridiagonalMatrixSecondDerivType2(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughSecondDerivativeType3:
			BuildingTridiagonalMatrixSecondDerivType3(x, n, h, matrix, d);
			break;
		case SplineRepresentationType::throughSecondDerivativeType4:
			BuildingTridiagonalMatrixSecondDerivType4(x, n, h, matrix, d);
			break;
	}

	if (bonCond == SplineRepresentationType::throughFirstDerivativeType3 || bonCond == SplineRepresentationType::throughSecondDerivativeType3) {
		alpha[1] = -matrix[1][2] / matrix[1][1];
		beta[1] = d[1] / matrix[1][1];

		for (int i = 2; i < n; i++) {
			alpha[i] = -matrix[i][i + 1] / (matrix[i][i] + matrix[i][i - 1] * alpha[i - 1]);
			beta[i] = (d[i] - matrix[i][i - 1] * beta[i - 1]) / (matrix[i][i] + matrix[i][i - 1] * alpha[i - 1]);
		}
	}
	else {
		alpha[0] = -matrix[0][1] / matrix[0][0];
		beta[0] = d[0] / matrix[0][0];

		for (int i = 1; i < n; i++) {
			alpha[i] = -matrix[i][i + 1] / (matrix[i][i] + matrix[i][i - 1] * alpha[i - 1]);
			beta[i] = (d[i] - matrix[i][i - 1] * beta[i - 1]) / (matrix[i][i] + matrix[i][i - 1] * alpha[i - 1]);
		}
	}

	for (int i = n - 1; i > 0; i--) {
		coefM[i] = alpha[i] * coefM[i + 1] + beta[i];
	}

	delete[] alpha;
	delete[] beta;
	delete[] d;
	for (int i = 0; i <= n; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}
#pragma endregion

#pragma region Building spline
/// <summary>
/// Построение сплайна от одной переменной
/// </summary>
double BuildingSpline(double* x, double* coefM, int n, double* h, double valueX, BuildingSplineType type) {
	int j = 0;

	for (int i = 0; i < n; i++) {
		if (valueX >= x[i] && valueX <= x[i + 1]) {
			j = i;
			break;
		}
	}

	double t = (valueX - x[j]) / h[j];

	if (type == BuildingSplineType::buildingSplineUsingFirstDerivative) {
		return ExactSolution(x[j]) * pow(1 - t, 2) * (1 + 2 * t) + ExactSolution(x[j + 1]) * pow(t, 2) * (3 - 2 * t) + coefM[j] * h[j] * t * pow(1 - t, 2) - coefM[j + 1] * h[j] * pow(t, 2) * (1 - t);
	}

	if (type == BuildingSplineType::buildingSplineUsingSecondDerivative) {
		return ExactSolution(x[j]) * (1 - t) + ExactSolution(x[j + 1]) * t - coefM[j] * t * (1 - t) * (2 - t) * pow(h[j], 2) / 6 + coefM[j + 1] * t * (pow(t, 2) - 1) * pow(h[j], 2) / 6;
	}
}
#pragma endregion