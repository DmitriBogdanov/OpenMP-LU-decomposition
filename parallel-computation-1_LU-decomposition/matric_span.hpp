#pragma once

#include "matrix.hpp"

// Save inverse of source span
// into dest span
template <typename T>
void span_inverse_LU(
	Matrix<T> &src, size_t src_i, size_t src_j, size_t rows, size_t cols,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	/// LET SIZE BE b = 2 SO THE INVERSE IS SIMPLE
	const auto a = T(1);// src(src_i, src_j);
	const auto b = T(0);// src(src_i, src_j + 1);
	const auto c = src(src_i + 1, src_j);
	const auto d = T(1);// src(src_i + 1, src_j + 1);

	const T inverseDet = T(1) / (a * d - b * c);

	dst(dst_i, dst_j) = inverseDet * d;
	dst(dst_i, dst_j + 1) = inverseDet * -b;
	dst(dst_i + 1, dst_j) = inverseDet * -c;
	dst(dst_i + 1, dst_j + 1) = inverseDet * a;
}

// Substract product of source1 and source2 spans
// from dest span
template <typename T>
void span_substract_product(
	Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	if (cols1 != rows2) throw std::runtime_error("span_substract_product(): incompatible matrices encountered.");

	for (size_t i = 0; i < rows1; ++i)
		for (size_t k = 0; k < cols1; ++k)
			for (size_t j = 0; j < cols2; ++j)
				dst(dst_i + i, dst_j + j) -= src1(src1_i + i, src1_j + k) * src2(src2_i + k, src2_j + j);
}

// Save product of source1 and source2 spans
// into dest span
template <typename T>
void span_set_product(
	Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	if (cols1 != rows2) {
		std::cout << "span_substract_product(): incompatible matrices encountered.";
		exit(1);
	}

	for (size_t i = 0; i < rows1; ++i)
		for (size_t j = 0; j < cols2; ++j)
			dst(dst_i + i, dst_j + j) = 0;

	for (size_t i = 0; i < rows1; ++i)
		for (size_t k = 0; k < cols1; ++k)
			for (size_t j = 0; j < cols2; ++j)
				dst(dst_i + i, dst_j + j) += src1(src1_i + i, src1_j + k) * src2(src2_i + k, src2_j + j);
}


template<typename T>
void inverse(
	Matrix<T>& src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	Matrix<T>& src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2)
{
	if (rows1 != rows2) {
		std::cout << "inverse_product(): matrix is not square.";
		exit(1);
	}

	for(size_t j = 0; j < cols2; ++j)
		for(size_t i = 0; i < rows2; ++i) 
			for(int k = i - 1; k >= 0; --k)
				src2(src2_i + i, src2_j + j) -= src2(src2_i + k, src2_j + j) * src1(src1_i + i, src1_j + k);
}
