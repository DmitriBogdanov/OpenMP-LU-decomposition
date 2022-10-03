#pragma once

#include "matric_parspan.hpp"


// # Parallel LU #
// - Sequential
// - No pivoting
// - Time complexity O(2/3 N^3)
template <typename T>
void LU_par(T *A, const size_t ROWS, const size_t COLS) {
	for (int i = 0; i < std::min(ROWS - 1, COLS); ++i) {
		// (1)
		const T inverseAii = T(1) / A[i * COLS + i];

		for (int j = i + 1; j < ROWS; ++j)
			A[j * COLS + i] *= inverseAii;

		// (2)
#pragma omp parallel for
		for (int j = i + 1; j < ROWS; ++j)
			for (int k = i + 1; k < COLS; ++k)
				A[j * COLS + k] -= A[j * COLS + i] * A[i * COLS + k];
	}
}


// # Parallel Block LU #
// - Sequential
// - No pivoting
// - Time complexity O(?)
template <typename T>
void LU_par_block(T *A, const size_t N, const size_t b) {
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