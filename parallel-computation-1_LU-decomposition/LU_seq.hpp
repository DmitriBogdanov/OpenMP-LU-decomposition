#pragma once

#include "matrix.hpp"

// LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(2/3 N^3)
template <typename T>
void LU_seq(const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U) {
	const size_t N = A.rows();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < N; ++j) {

			if (j < i)
				L(j, i) = 0;
			else {
				L(j, i) = A(j, i);

				for (size_t k = 0; k < i; ++k) {
					L(j, i) = L(j, i) - L(j, k) * U(k, i);
				}
			}
		}

		for (size_t j = 0; j < N; j++) {
			if (j < i)
				U(i, j) = 0;
			else if (j == i)
				U(i, j) = 1;
			else {
				U(i, j) = A(i, j) / L(i, i);
				for (size_t k = 0; k < i; ++k) {
					U(i, j) = U(i, j) - ((L(i, k) * U(k, j)) / L(i, i));
				}
			}
		}
	}
}