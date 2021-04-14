#include "InitialFunctions.h"

#pragma region Enumerations
/// <summary>
/// Производная, через которую представляем сплайн, и тип краевых условий
/// </summary>
enum class SplineRepresentationType {
	throughFirstDerivativeType1,
	throughFirstDerivativeType2,
	throughFirstDerivativeType3,
	throughFirstDerivativeType4,
	throughSecondDerivativeType1,
	throughSecondDerivativeType2,
	throughSecondDerivativeType3,
	throughSecondDerivativeType4
};

/// <summary>
/// Тип производной, который используется для построение сплайна 
/// </summary>
enum class BuildingSplineType {
	buildingSplineUsingFirstDerivative,
	buildingSplineUsingSecondDerivative
};

/// <summary>
/// цикл и фиксируемая переменная 
/// </summary>
enum class FixedVariableType {
	firstCycleFixedVariableX,
	firstCycleFixedVariableY,
	secondCycleFixedVariableY
};
#pragma endregion

#pragma region Support function
/// <summary>
/// Вычисление функции phi
/// </summary>
void CalcPhi(double* phi, double h, double tau, BuildingSplineType type) {
	if (type == BuildingSplineType::buildingSplineUsingFirstDerivative) {
		phi[0] = pow(1 - tau, 2) * (1 + 2 * tau);
		phi[1] = pow(tau, 2) * (3 - 2 * tau);
		phi[2] = tau * pow(1 - tau, 2) * h;
		phi[3] = pow(tau, 2) * (tau - 1) * h;
	}

	if (type == BuildingSplineType::buildingSplineUsingSecondDerivative) {
		phi[0] = 1 - tau;
		phi[1] = tau;
		phi[2] = tau * (1 - tau) * (2 - tau) * h * h;
		phi[3] = tau * (pow(tau, 2) - 1) * h * h;
	}
}

/// <summary>
/// Вычисление произведения матриц
/// </summary>
double MatrixMultiplication(double* matrixPhiT, double** matrixF, double* matrixPhiU) {
	double* temp = new double[4];
	double result = 0;

	for (int i = 0; i < 4; i++) {
		temp[i] = 0;
		for (int j = 0; j < 4; j++) {
			temp[i] += matrixPhiT[j] * matrixF[j][i];
		}
	}

	for (int i = 0; i < 4; i++) {
		result += temp[i] * matrixPhiU[i];
	}

	delete[] temp;

	return result;
}
#pragma endregion

#pragma region Boundary conditions
/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа I
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType1(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
	matrix[0][0] = 1;												// coef for m_0
	matrix[0][1] = 0;												// coef for m_1
	matrix[n][n - 1] = 0;											// coef for m_N-1
	matrix[n][n] = 1;												// coef for m_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i] / (h[i - 1] + h[i]);				//coef for m_i-1 -  \lambda_i
		matrix[i][i] = 2;											//coef for m_i
		matrix[i][i + 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for m_i+1 -  \mu_i
	}

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[0] = PartialDerivative01(fixedVariableValue, x[0]);
		d[n] = PartialDerivative01(fixedVariableValue, x[n]);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / (h[i - 1] + h[i]) / h[i] + h[i] 
				* (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[0] = PartialDerivative10(x[0], fixedVariableValue);
		d[n] = PartialDerivative10(x[n], fixedVariableValue);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i] + h[i]
				* (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}

	if(fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[0] = MixedPartialDerivative11(x[0], fixedVariableValue);
		d[n] = MixedPartialDerivative11(x[n], fixedVariableValue);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (coefM[i + 1] - coefM[i]) / (h[i - 1] + h[i]) / h[i] + h[i] * (coefM[i] - coefM[i - 1]) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа II
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType2(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
	matrix[0][0] = 2;												// coef for m_0
	matrix[0][1] = 1;												// coef for m_1
	matrix[n][n - 1] = 1;											// coef for m_N-1
	matrix[n][n] = 2;												// coef for m_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i] / (h[i - 1] + h[i]);				//coef for m_i-1 -  \lambda_i
		matrix[i][i] = 2;											//coef for m_i
		matrix[i][i + 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for m_i+1 -  \mu_i
	}

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[0] = 3 * (ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / h[0] - 0.5 * h[0] * PartialDerivative02(fixedVariableValue, x[0]);
		d[n] = 3 * (ExactSolution(fixedVariableValue, x[n]) - ExactSolution(fixedVariableValue, x[n - 1])) / h[n - 1] - 0.5 * h[n - 1] * PartialDerivative02(fixedVariableValue, x[n]);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / (h[i - 1] + h[i]) / h[i] + h[i] 
				* (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[0] = 3 * (ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / h[0] - 0.5 * h[0] * PartialDerivative20(x[0], fixedVariableValue);
		d[n] = 3 * (ExactSolution(x[n], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / h[n - 1] - 0.5 * h[n - 1] * PartialDerivative20(x[n], fixedVariableValue);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i] + h[i] 
				* (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}

	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[0] = 3 * (coefM[1] - coefM[0]) / h[0] - 0.5 * h[0] * MixedPartialDerivative22(x[0], fixedVariableValue);
		d[n] = 3 * (coefM[n] - coefM[n - 1]) / h[n - 1] - 0.5 * h[n - 1] * MixedPartialDerivative22(x[n], fixedVariableValue);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (coefM[i + 1] - coefM[i]) / (h[i - 1] + h[i]) / h[i] + h[i]	* (coefM[i] - coefM[i - 1]) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа III
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType3(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
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

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[n] = 3 * (h[n - 1] * (ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / (h[n - 1] + h[0]) / h[0] + h[0]
			* (ExactSolution(fixedVariableValue, x[0]) - ExactSolution(fixedVariableValue, x[n - 1])) / (h[n - 1] + h[0]) / h[n - 1]);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / (h[i - 1] + h[i]) / h[i] + h[i]
				* (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[n] = 3 * (h[n - 1] * (ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / (h[n - 1] + h[0]) / h[0] + h[0]
			* (ExactSolution(x[0], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / (h[n - 1] + h[0]) / h[n - 1]);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i] + h[i]
				* (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[n] = 3 * (h[n - 1] * (coefM[1] - coefM[0]) / (h[n - 1] + h[0]) / h[0] + h[0] * (coefM[0] - coefM[n - 1]) / (h[n - 1] + h[0]) / h[n - 1]);

		for (int i = 1; i < n; i++) {
			d[i] = 3 * (h[i - 1] * (coefM[i + 1] - coefM[i]) / (h[i - 1] + h[i]) / h[i] + h[i] * (coefM[i] - coefM[i - 1]) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через первую производную, с краевыми условиями типа IV
/// </summary>
void BuildingTridiagonalMatrixFirstDerivType4(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
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

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[0] = 2 * ((ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / h[0] - pow(h[0] / h[1], 2) 
			* (ExactSolution(fixedVariableValue, x[2]) - ExactSolution(fixedVariableValue, x[1])) / h[1]);
		d[1] = (h[0] / (h[0] + h[1]) + 2 * h[0] / h[1]) * (ExactSolution(fixedVariableValue, x[2]) - ExactSolution(fixedVariableValue, x[1])) / h[1] + h[1]
			* (ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / (h[0] + h[1]) / h[0];
		d[n - 1] = h[n - 2] * (ExactSolution(fixedVariableValue, x[n]) - ExactSolution(fixedVariableValue, x[n - 1])) / (h[n - 2] + h[n - 1]) / h[n - 1] + (h[n - 1] / (h[n - 2] + h[n - 1]) 
			+ 2 * h[n - 1] / h[n - 2]) * (ExactSolution(fixedVariableValue, x[n - 1]) - ExactSolution(fixedVariableValue, x[n - 2])) / h[n - 2];
		d[n] = 2 * ((ExactSolution(fixedVariableValue, x[n]) - ExactSolution(fixedVariableValue, x[n - 1])) / h[n - 1] - pow(h[n - 1] / h[n - 2], 2) 
			* (ExactSolution(fixedVariableValue, x[n - 1]) - ExactSolution(fixedVariableValue, x[n - 2])) / h[n - 2]);

		for (int i = 2; i < n - 1; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / (h[i - 1] + h[i]) / h[i] + h[i]
				* (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[0] = 2 * ((ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / h[0] - pow(h[0] / h[1], 2) 
			* (ExactSolution(x[2], fixedVariableValue) - ExactSolution(x[1], fixedVariableValue)) / h[1]);
		d[1] = (h[0] / (h[0] + h[1]) + 2 * h[0] / h[1]) * (ExactSolution(x[2], fixedVariableValue) - ExactSolution(x[1], fixedVariableValue)) / h[1] + h[1] 
			* (ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / (h[0] + h[1]) / h[0];
		d[n - 1] = h[n - 2] * (ExactSolution(x[n], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / (h[n - 2] + h[n - 1]) / h[n - 1] + (h[n - 1] / (h[n - 2] + h[n - 1]) 
			+ 2 * h[n - 1] / h[n - 2]) * (ExactSolution(x[n - 1], fixedVariableValue) - ExactSolution(x[n - 2], fixedVariableValue)) / h[n - 2];
		d[n] = 2 * ((ExactSolution(x[n], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / h[n - 1] - pow(h[n - 1] / h[n - 2], 2) 
			* (ExactSolution(x[n - 1], fixedVariableValue) - ExactSolution(x[n - 2], fixedVariableValue)) / h[n - 2]);

		for (int i = 2; i < n - 1; i++) {
			d[i] = 3 * (h[i - 1] * (ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i] + h[i] 
				* (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[0] = 2 * ((coefM[1] - coefM[0]) / h[0] - pow(h[0] / h[1], 2) * (coefM[2] - coefM[1]) / h[1]);
		d[1] = (h[0] / (h[0] + h[1]) + 2 * h[0] / h[1]) * (coefM[2] - coefM[1]) / h[1] + h[1] * (coefM[1] - coefM[0]) / (h[0] + h[1]) / h[0];
		d[n - 1] = h[n - 2] * (coefM[n] - coefM[n - 1]) / (h[n - 2] + h[n - 1]) / h[n - 1] + (h[n - 1] / (h[n - 2] + h[n - 1]) + 2 * h[n - 1] / h[n - 2]) * (coefM[n - 1] - coefM[n - 2]) / h[n - 2];
		d[n] = 2 * ((coefM[n] - coefM[n - 1]) / h[n - 1] - pow(h[n - 1] / h[n - 2], 2)	* (coefM[n - 1] - coefM[n - 2]) / h[n - 2]);

		for (int i = 2; i < n - 1; i++) {
			d[i] = 3 * (h[i - 1] * (coefM[i + 1] - coefM[i]) / (h[i - 1] + h[i]) / h[i] + h[i] * (coefM[i] - coefM[i - 1]) / (h[i - 1] + h[i]) / h[i - 1]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа I
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType1(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
	matrix[0][0] = 2;												// coef for M_0
	matrix[0][1] = 1;												// coef for M_1
	matrix[n][n - 1] = 1;											// coef for M_N-1
	matrix[n][n] = 2;												// coef for M_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for M_i-1 - \mu_i
		matrix[i][i] = 2;											//coef for M_i
		matrix[i][i + 1] = h[i] / (h[i - 1] + h[i]);				//coef for M_i+1 -  \lambda_i
	}

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[0] = 6 * ((ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / h[0] - PartialDerivative01(fixedVariableValue, x[0])) / h[0];
		d[n] = 6 * (PartialDerivative01(fixedVariableValue, x[n]) - (ExactSolution(fixedVariableValue, x[n]) - ExactSolution(fixedVariableValue, x[n - 1])) / h[n - 1]) / h[n - 1];

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / h[i] 
				- (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[0] = 6 * ((ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / h[0] - PartialDerivative10(x[0], fixedVariableValue)) / h[0];
		d[n] = 6 * (PartialDerivative10(x[n], fixedVariableValue) - (ExactSolution(x[n], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / h[n - 1]) / h[n - 1];

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / h[i] 
				- (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[0] = 6 * ((coefM[1] - coefM[0]) / h[0] - MixedPartialDerivative11(x[0], fixedVariableValue)) / h[0];
		d[n] = 6 * (MixedPartialDerivative11(x[n], fixedVariableValue) - (coefM[n] - coefM[n - 1]) / h[n - 1]) / h[n - 1];

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((coefM[i + 1] - coefM[i]) / h[i] - (coefM[i] - coefM[i - 1]) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа II
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType2(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
	matrix[0][0] = 1;												// coef for M_0
	matrix[0][1] = 0;												// coef for M_1
	matrix[n][n - 1] = 0;											// coef for M_N-1
	matrix[n][n] = 1;												// coef for M_N

	for (int i = 1; i < n; i++) {
		matrix[i][i - 1] = h[i - 1] / (h[i - 1] + h[i]);			//coef for M_i-1 - \mu_i
		matrix[i][i] = 2;											//coef for M_i
		matrix[i][i + 1] = h[i] / (h[i - 1] + h[i]);				//coef for M_i+1 -  \lambda_i
	}

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[0] = PartialDerivative02(fixedVariableValue, x[0]);
		d[n] = PartialDerivative02(fixedVariableValue, x[n]);

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / h[i] 
				- (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[0] = PartialDerivative20(x[0], fixedVariableValue);
		d[n] = PartialDerivative20(x[n], fixedVariableValue);

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / h[i] 
				- (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[0] = MixedPartialDerivative22(x[0], fixedVariableValue);
		d[n] = MixedPartialDerivative22(x[n], fixedVariableValue);

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((coefM[i + 1] - coefM[i]) / h[i] - (coefM[i] - coefM[i - 1]) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа III
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType3(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
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

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[n] = 6 * ((ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / h[0] 
			- (ExactSolution(fixedVariableValue, x[0]) - ExactSolution(fixedVariableValue, x[n - 1])) / h[n - 1]) / (h[n - 1] + h[0]);

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / h[i] 
				- (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[n] = 6 * ((ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / h[0] 
			- (ExactSolution(x[0], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / h[n - 1]) / (h[n - 1] + h[0]);

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / h[i] 
				- (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[n] = 6 * ((coefM[1] - coefM[0]) / h[0] - (coefM[0] - coefM[n - 1]) / h[n - 1]) / (h[n - 1] + h[0]);

		for (int i = 1; i < n; i++) {
			d[i] = 6 * ((coefM[i + 1] - coefM[i]) / h[i] - (coefM[i] - coefM[i - 1]) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
}

/// <summary>
/// Построение диагональной матрицы сплайна, представленного через вторую производную, с краевыми условиями типа IV
/// </summary>
void BuildingTridiagonalMatrixSecondDerivType4(double* x, int n, double* h, double** matrix, double* d, double fixedVariableValue, FixedVariableType fixedVariable, double* coefM) {
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

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		d[1] = 6 * h[1] * ((ExactSolution(fixedVariableValue, x[2]) - ExactSolution(fixedVariableValue, x[1])) / h[1] 
			- (ExactSolution(fixedVariableValue, x[1]) - ExactSolution(fixedVariableValue, x[0])) / h[0]) / pow(h[0] + h[1], 2);
		d[n - 1] = 6 * h[n - 2] * ((ExactSolution(fixedVariableValue, x[n]) - ExactSolution(fixedVariableValue, x[n - 1])) / h[n - 1] 
			- (ExactSolution(fixedVariableValue, x[n - 1]) - ExactSolution(fixedVariableValue, x[n - 2])) / h[n - 2]) / pow(h[n - 2] + h[n - 1], 2);

		for (int i = 2; i < n - 1; i++) {
			d[i] = 6 * ((ExactSolution(fixedVariableValue, x[i + 1]) - ExactSolution(fixedVariableValue, x[i])) / h[i] 
				- (ExactSolution(fixedVariableValue, x[i]) - ExactSolution(fixedVariableValue, x[i - 1])) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		d[1] = 6 * h[1] * ((ExactSolution(x[2], fixedVariableValue) - ExactSolution(x[1], fixedVariableValue)) / h[1] 
			- (ExactSolution(x[1], fixedVariableValue) - ExactSolution(x[0], fixedVariableValue)) / h[0]) / pow(h[0] + h[1], 2);
		d[n - 1] = 6 * h[n - 2] * ((ExactSolution(x[n], fixedVariableValue) - ExactSolution(x[n - 1], fixedVariableValue)) / h[n - 1] 
			- (ExactSolution(x[n - 1], fixedVariableValue) - ExactSolution(x[n - 2], fixedVariableValue)) / h[n - 2]) / pow(h[n - 2] + h[n - 1], 2);

		for (int i = 2; i < n - 1; i++) {
			d[i] = 6 * ((ExactSolution(x[i + 1], fixedVariableValue) - ExactSolution(x[i], fixedVariableValue)) / h[i] 
				- (ExactSolution(x[i], fixedVariableValue) - ExactSolution(x[i - 1], fixedVariableValue)) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		d[1] = 6 * h[1] * ((coefM[2] - coefM[1]) / h[1]	- (coefM[1] - coefM[0]) / h[0]) / pow(h[0] + h[1], 2);
		d[n - 1] = 6 * h[n - 2] * ((coefM[n] - coefM[n - 1]) / h[n - 1]	- (coefM[n - 1] - coefM[n - 2]) / h[n - 2]) / pow(h[n - 2] + h[n - 1], 2);

		for (int i = 2; i < n - 1; i++) {
			d[i] = 6 * ((coefM[i + 1] - coefM[i]) / h[i] - (coefM[i] - coefM[i - 1]) / h[i - 1]) / (h[i - 1] + h[i]);
		}
	}
}
#pragma endregion

#pragma region Sweep method
/// <summary>
/// Метод прогонки
/// </summary>
void SweepMethod(double* x, double* coefM, int n, double* h, SplineRepresentationType bonCond, double fixedVariableValue, FixedVariableType fixedVariable, double* tempCoefM) {
	double* alpha = new double[n]();
	double* beta = new double[n]();

	double* d = new double[n + 1]();
	double** matrix = new double* [n + 1];
	for (int i = 0; i <= n; i++) {
		matrix[i] = new double[n + 1]();
	}

	if (fixedVariable == FixedVariableType::firstCycleFixedVariableX) {
		coefM[0] = ExactSolution(fixedVariableValue, x[0]);
		coefM[n] = ExactSolution(fixedVariableValue, x[n]);
	}
	
	if (fixedVariable == FixedVariableType::firstCycleFixedVariableY) {
		coefM[0] = ExactSolution(x[0], fixedVariableValue);
		coefM[n] = ExactSolution(x[n], fixedVariableValue);
	}
	
	if (fixedVariable == FixedVariableType::secondCycleFixedVariableY) {
		coefM[0] = tempCoefM[0];
		coefM[n] = tempCoefM[n];
	}
	
	switch (bonCond) {
		case SplineRepresentationType::throughFirstDerivativeType1:
			BuildingTridiagonalMatrixFirstDerivType1(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughFirstDerivativeType2:
			BuildingTridiagonalMatrixFirstDerivType2(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughFirstDerivativeType3:
			BuildingTridiagonalMatrixFirstDerivType3(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughFirstDerivativeType4:
			BuildingTridiagonalMatrixFirstDerivType4(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughSecondDerivativeType1:
			BuildingTridiagonalMatrixSecondDerivType1(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughSecondDerivativeType2:
			BuildingTridiagonalMatrixSecondDerivType2(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughSecondDerivativeType3:
			BuildingTridiagonalMatrixSecondDerivType3(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
			break;
		case SplineRepresentationType::throughSecondDerivativeType4:
			BuildingTridiagonalMatrixSecondDerivType4(x, n, h, matrix, d, fixedVariableValue, fixedVariable, tempCoefM);
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

#pragma region BuildingSpline
/// <summary>
/// Построение сплайна от двух переменных
/// </summary>
double BuildingSpline(double* x, double* y, double** coefM10, double** coefM01, double** coefM11, int nX, int nY, double* hX, double* hY, double valueX, double valueY,	BuildingSplineType type) {
	int k = 0, l = 0;

	for (int i = 0; i < nX; i++) {
		if (valueX >= x[i] && valueX <= x[i + 1]) {
			k = i;
			break;
		}
	}

	for (int j = 0; j < nY; j++) {
		if (valueY >= y[j] && valueY <= y[j + 1]) {
			l = j;
			break;
		}
	}

	double t = (valueX - x[k]) / hX[k];
	double u = (valueY - y[l]) / hY[l];

	double** matrixF = new double* [4];
	for (int i = 0; i < 4; i++) {
		matrixF[i] = new double[4];
	}

	double* matrixPhiT = new double[4];
	double* matrixPhiU = new double[4];
	double result;
	int val = 1;

	if (type == BuildingSplineType::buildingSplineUsingSecondDerivative) {
		val = 6;
	}

	matrixF[0][0] = ExactSolution(x[k],y[l]);
	matrixF[0][1] = ExactSolution(x[k + 1],y[l]);
	matrixF[0][2] = coefM10[k][l] / val;
	matrixF[0][3] = coefM10[k + 1][l] / val;

	matrixF[1][0] = ExactSolution(x[k], y[l + 1]);
	matrixF[1][1] = ExactSolution(x[k + 1], y[l + 1]);
	matrixF[1][2] = coefM10[k][l + 1] / val;
	matrixF[1][3] = coefM10[k + 1][l + 1] / val;

	matrixF[2][0] = coefM01[k][l] / val;
	matrixF[2][1] = coefM01[k + 1][l] / val;
	matrixF[2][2] = coefM11[k][l] / val / val;
	matrixF[2][3] = coefM11[k + 1][l] / val / val;

	matrixF[3][0] = coefM01[k][l + 1] / val;
	matrixF[3][1] = coefM01[k + 1][l + 1] / val;
	matrixF[3][2] = coefM11[k][l + 1] / val / val;
	matrixF[3][3] = coefM11[k + 1][l + 1] / val / val;

	CalcPhi(matrixPhiT, hX[k], t, type);
	CalcPhi(matrixPhiU, hY[l], u, type);
	result = MatrixMultiplication(matrixPhiT, matrixF, matrixPhiU);

	for (int i = 0; i < 4; i++) {
		delete[] matrixF[i];
	}
	delete[] matrixF;
	delete[] matrixPhiT;
	delete[] matrixPhiU;

	return result;
}
#pragma endregion