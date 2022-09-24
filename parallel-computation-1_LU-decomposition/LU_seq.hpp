#pragma once

#include "matrix.hpp"
#include "matric_span.hpp"


// LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(2/3 N^3)
template <typename T>
void LU_seq(Matrix<T> &A) {
	const size_t ROWS = A.rows();
	const size_t COLS = A.cols();

	for (size_t i = 0; i < std::min(ROWS - 1, COLS); ++i) {
		// (1)
		const T inverseAii = T(1) / A(i, i);

		for (size_t j = i + 1; j < ROWS; ++j)
			A(j, i) *= inverseAii;

		// (2)
		for (size_t j = i + 1; j < ROWS; ++j)
			for (size_t k = i + 1; k < COLS; ++k)
				A(j, k) -= A(j, i) * A(i, k);
	}
}

// Block LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(?)
template <typename T>
void LU_seq_block(Matrix<T>& A, const size_t b) {
	// ROWS = COLS like in the literature
	const size_t N = A.rows();

	auto inverseL22 = Matrix<T>(b, b); // temp storage for inversion result
	auto temp = Matrix<T>(b, N - b);

	for (size_t i = 0; i < N; i += b) {
		// (1)
		// Find LU decomposition of block
		const size_t rows = N - i;
		const size_t cols = b;

		for (size_t _i = 0; _i < std::min(rows - 1, cols); ++_i) {
			const T inverseAii = T(1) / A(_i + i, _i + i);

			for (size_t _j = _i + 1; _j < rows; ++_j)
				A(_j + i, _i + i) *= inverseAii;

			for (size_t _j = _i + 1; _j < rows; ++_j)
				for (size_t _k = _i + 1; _k < cols; ++_k)
					A(_j + i, _k + i) -= A(_j + i, _i + i) * A(_i + i, _k + i);
		}

		// (2)
		// Solve systems L22 * U23 = A23
		// to get U23 = L22^-1 * A23
		span_set_product_inverseL_by_self(
			A,
			i, i,
			b, b,
			A,
			i, i + b,
			b, N - b - i
		);
	
		
		// (3)
		// A33 -= A32 * A23
		span_substract_product(
			A,
			i + b, i,
			N - b - i, b,
			A,
			i, i + b,
			b, N - b - i,
			A,
			i + b,
			i + b
		);
	}
}

// Splits LU saved as a single matrix into separate objects
template <typename T>
void split_into_LU(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U) {
	const auto MIN = std::min(A.rows(), A.cols());

	// L (ROWS x ROWS)
	for (size_t i = 0; i < A.rows(); ++i) {
		for (size_t j = 0; j < std::min(i, MIN); ++j) L(i, j) = A(i, j);
		if (i < MIN) L(i, i) = 1;
		for (size_t j = i + 1; j < MIN; ++j) L(i, j) = 0;
	}

	// U (ROWS x COLS)
	for (size_t i = 0; i < MIN; ++i) {
		for (size_t j = 0; j < i; ++j) U(i, j) = 0;
		for (size_t j = i; j < A.cols(); ++j) U(i, j) = A(i, j);
	}
}