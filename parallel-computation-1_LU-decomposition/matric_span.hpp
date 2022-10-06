#pragma once

#include "matrix.hpp"


// Copy a block into another
template<typename T>
inline void span_copy(
	Matrix<T> &src, size_t src_i, size_t src_j, size_t rows, size_t cols,
	Matrix<T> &dst, size_t dst_i, size_t dst_j, bool transposed = false) {

	// Copy each row with memcpy()
	for (size_t i = 0; i < rows; ++i)
		memcpy(&dst(dst_i + i, dst_j), &src(src_i + i, src_j), cols * sizeof(T));
}


// Copy a block into another
// when 'src' and 'dst' matrices have the same width we can use a single memcpy() to
// copy multiple rows
template<typename T>
inline void span_copy_same_cols(
	const Matrix<T> &src, size_t src_i, size_t src_j, size_t rows, size_t cols,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {
	
	assert(src.cols() == dst.cols() && "span_copy_same_cols(): incompatible matrices encountered.");

	// Copy all rows with a single memcpy()
	memcpy(&dst(dst_i, dst_j), &src(src_i, src_j), rows * cols * sizeof(T));
}


// LU decompose src
template<typename T>
inline void span_LU_decomposition(
	Matrix<T> &src, size_t src_i, size_t src_j, size_t rows, size_t cols) {

	for (size_t _i = 0; _i < std::min(rows - 1, cols); ++_i) {
		const T inverseAii = T(1) / src(src_i + _i, src_j + _i);

		for (size_t _j = _i + 1; _j < rows; ++_j)
			src(src_i + _j, src_j + _i) *= inverseAii;

		for (size_t _j = _i + 1; _j < rows; ++_j)
			for (size_t _k = _i + 1; _k < cols; ++_k)
				src(src_i + _j, src_j + _k) -= src(src_i + _j, src_j + _i) * src(src_i + _i, src_j + _k);
	}
}


// Solve
// into dest span
template<typename T>
inline void span_set_product_inverseL_by_self(
	const Matrix<T>& src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	Matrix<T>& src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2)
{
	assert(rows1 == rows2 && "inverse_product(): incompatible matrices encountered.");

	for (size_t j = 0; j < cols2; ++j)
		for (size_t i = 0; i < rows2; ++i)
			for (int k = i - 1; k >= 0; --k)
				src2(src2_i + i, src2_j + j) -= src2(src2_i + k, src2_j + j) * src1(src1_i + i, src1_j + k);
}

// Substract product of source1 and source2 spans
// from dest span
template <typename T>
inline void span_substract_product(
	const Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	const Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	assert(cols1 == rows2 && "span_substract_product(): incompatible matrices encountered.");

	for (size_t i = 0; i < rows1; ++i)
		for (size_t k = 0; k < cols1; ++k)
			for (size_t j = 0; j < cols2; ++j)
				dst(dst_i + i, dst_j + j) -= src1(src1_i + i, src1_j + k) * src2(src2_i + k, src2_j + j);
}


template <typename T>
inline void span_substract_product_block(
	const Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	const Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	assert(cols1 == rows2 && "span_substract_product(): incompatible matrices encountered.");

	matrix_mult_wiki_block(&(src1._data[0]), &(src2._data[0]), &(dst._data[0]), rows1, cols2, cols1);
}
