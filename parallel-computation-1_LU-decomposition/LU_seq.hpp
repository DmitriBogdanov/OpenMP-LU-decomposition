#pragma once

#include <memory>

#include "matrix.hpp"
#include "matric_span.hpp"


// LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(2/3 N^3)
template <typename T>
void LU_seq(T *A, const size_t ROWS, const size_t COLS) {
	for (size_t i = 0; i < std::min(ROWS - 1, COLS); ++i) {
		// (1)
		const T inverseAii = T(1) / A[i * COLS + i];

		for (size_t j = i + 1; j < ROWS; ++j)
			A[j * COLS + i] *= inverseAii;

		// (2)
		for (size_t j = i + 1; j < ROWS; ++j)
			for (size_t k = i + 1; k < COLS; ++k)
				A[j * COLS + k] -= A[j * COLS + i] * A[i * COLS + k];
	}
}

// Block LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(?)
template <typename T>
void LU_seq_block(T *A, const size_t N, const size_t b) {
	const size_t total_length = N * b + b * (N - b);

	T* const buffer = new T[total_length * sizeof(T)];

	T *A_22;
	T *A_32;
	T *A_23; // NOTE: A_23 is col-major!

	for (size_t i = 0; i < N; i += b) {
		// Adjust pointers. Blocks A_22 -> A_32 -> A_23 are stored in corresponding order
		const size_t rows_22 = b;
		const size_t cols_22 = b;

		const size_t rows_32 = N - b - i;
		const size_t cols_32 = b;

		const size_t rows_23 = b;
		const size_t cols_23 = N - b - i;

		A_22 = buffer;
		A_32 = A_22 + rows_22 * cols_22;
		A_23 = A_32 + rows_32 * cols_32;

		// (1)
		// Find LU decomposition of block (A22 & A32)
		span_copy_rm_to_rm(
			// source
			A, N, N,
			i, i,
			rows_22 + rows_32, cols_22,
			// dest
			A_22, rows_22 + rows_32, cols_22,
			0, 0
		);

		LU_seq(A_22, rows_22 + rows_32, cols_22);

		

		span_copy_rm_to_rm(
			// source
			A_22, rows_22 + rows_32, cols_22,
			0, 0,
			rows_22 + rows_32, cols_22,
			// dest
			A, N, N,
			i, i
		);

		// (2)
		// Solve (N - b - i) systems L22 * x = A23
		// to get A23 = L22^-1 * A23
		span_copy_rm_to_cm(
			// source
			A, N, N,
			i, i + b,
			b, N - b - i,
			// dest
			A_23, rows_23, cols_23,
			0, 0
		);

		block_get_U23(
			A_22, rows_22, cols_22,
			A_23, rows_23, cols_23
		);

		span_copy_cm_to_rm(
			// source
			A_23, rows_23, cols_23,
			0, 0,
			rows_23, cols_23,
			// dest
			A, N, N,
			i, i + b
		);

		// (3)
		// A33 -= A32 * A23
		block_substract_product(
			// source 1
			A_32, rows_32, cols_32,
			// source 2
			A_23, rows_23, cols_23,
			// dest
			A, N, N,
			i + b, i + b
		);
	}

	delete[] buffer;
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

template <typename T>
Matrix<T> verify_LU(const Matrix<T> &A, const Matrix<T> &initial_matrix) {
	Matrix<T> res(initial_matrix);

	for (size_t i = 0; i < A.rows(); ++i)
		for (size_t j = 0; j < A.cols(); ++j)
			for (size_t k = 0; k <= std::min(j, i); ++k)
				if (k == i)
					res(i, j) -= A(i, j);
				else
					res(i, j) -= A(i, k) * A(k, j);

	return res;
}