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
		// As a result we U23 = L22^-1 * A23
		span_get_U23(
			A,
			i, i,
			b, b,
			A,
			i, i + b,
			b, N - b - i
		);
	
		
		// (3)
		// <block> -= <block 1> * <block 2>
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
	// L (ROWS x ROWS)
	for (size_t i = 0; i < A.rows(); ++i) {
		for (size_t j = 0; j < i; ++j) L(i, j) = A(i, j);
		L(i, i) = 1;
		for (size_t j = i + 1; j < A.rows(); ++j) L(i, j) = 0;
	}

	// U (ROWS x COLS)
	for (size_t i = 0; i < A.rows(); ++i) {
		for (size_t j = 0; j < i; ++j) U(i, j) = 0;
		for (size_t j = i; j < A.cols(); ++j) U(i, j) = A(i, j);
	}
}


// The only difference bethween method 2_9 and method 2_4 is 'if (i < COLS)'
// which seemingly doesn't have any effect on the behaviour since the condition is always true
/*
template <typename T>
void LU_seq_2_9(Matrix<T>& A) {
	const size_t ROWS = A.rows();
	const size_t COLS = A.cols();

	constexpr BLOCK_SIZE = 2;

	for (size_t i = 0; i < std::min(ROWS - 1, COLS); i += BLOCK_SIZE) {
		const T inverseAii = 1. / A(i, i);

		for (size_t j = i + 1; j < ROWS; ++j)
			A(j, i) *= inverseAii;

		if (i < COLS)
			for (size_t j = i + 1; j < ROWS; ++j)
				for (size_t k = i + 1; k < COLS; ++k)
					A(j, k) -= A(j, i) * A(i, k);
	}
}*/


//// (2)
//		// Save L22^-1
//span_inverse_LU(
//	A,
//	i, i,
//	b, b,
//	inverseL22,
//	0, 0
//);
//
//// Copy a part of 'A' into temp since source and dest of matrix multiplication can't overlap
//for (size_t _i = 0; _i < b; ++_i)
//	for (size_t _j = 0; _j < N - i - b; ++_j)
//		temp(_i, _j) = A(i + _i, i + b + _j);
//
/////std::cout << " i = " << i << "\nA = \n" << A << "\ntemp = \n" << temp << "\n\n";
//
//// <block> = L22^-1 * <block>
//span_set_product(
//	inverseL22,
//	0, 0,
//	b, b,
//	temp,
//	0, 0,
//	b, N - i - b,
//	A,
//	i, i + b
//);