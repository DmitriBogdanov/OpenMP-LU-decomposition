#pragma once

#include "matrix.hpp"
#include "static_timer.hpp"

// LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(2/3 N^3)
template <typename T>
void LU_seq(const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U) {
	const size_t N = A.rows();
	StaticTimer::start();
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
	std::cout << "LU decomposition time: " << StaticTimer::end() / 1000 << "(sec) \n";
}

// Block LU decomposition
// - Sequential
// - No pivoting
// - Time complexity O(?)
template <typename T>
void LU_seq_block(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U) {
	const size_t I = A.rows(), J = A.cols();
	Matrix<T> LU(A);

	size_t external_iters = (I - 1 > J) ? J : I - 1 ;

	T divideAii;
	StaticTimer::start();
	for (size_t i = 0; i < external_iters; ++i) {
		divideAii = 1 / LU(i, i);
		for (size_t j = i + 1; j < I; ++j) {
			LU(j, i) = LU(j, i) * divideAii;
		}
		if (i < J)
			for (size_t j = i + 1; j < I; ++j)
				for (size_t k = i + 1; k < J; ++k)
					LU(j, k) = LU(j, k) - LU(j, i) * LU(i, k);
	}
	std::cout << "Block LU decomposition time: " << StaticTimer::end() / 1000 << "(sec) \n";

	for (size_t i = 0; i < I; ++i) {
		L(i, i) = 1;
		U(i, i) = LU(i, i);
		for (size_t j = i + 1; j < J; ++j) {
			L(j, i) = LU(j, i);
			U(i, j) = LU(i, j);
		};
	}
}
