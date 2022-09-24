#pragma once

#include "matrix.hpp"

// Solve 
// into dest span
template<typename T>
void span_set_product_inverseL_by_self(
	const Matrix<T>& src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	Matrix<T>& src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2)
{
	if (rows1 != rows2) {
		std::cout << "inverse_product(): incompatible matrices encountered.";
		exit(1);
	}

	for (size_t j = 0; j < cols2; ++j)
		for (size_t i = 0; i < rows2; ++i)
			for (int k = i - 1; k >= 0; --k)
				src2(src2_i + i, src2_j + j) -= src2(src2_i + k, src2_j + j) * src1(src1_i + i, src1_j + k);
}

// Substract product of source1 and source2 spans
// from dest span
template <typename T>
void span_substract_product(
	const Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	const Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	if (cols1 != rows2) throw std::runtime_error("span_substract_product(): incompatible matrices encountered.");

	for (size_t i = 0; i < rows1; ++i)
		for (size_t k = 0; k < cols1; ++k)
			for (size_t j = 0; j < cols2; ++j)
				dst(dst_i + i, dst_j + j) -= src1(src1_i + i, src1_j + k) * src2(src2_i + k, src2_j + j);
}