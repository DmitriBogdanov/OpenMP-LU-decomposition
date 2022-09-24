#pragma once

#include "matrix.hpp"

// Solve 
// into dest span
template<typename T>
void span_set_product_inverseL_by_self(Matrix<T>& A, const Matrix<T>& L, Matrix<T>& U, size_t external_i)
{
	//static_assert(L.rows() != U.rows(), "inverse_product(): incompatible matrices encountered.");


	for (size_t j = 0; j < U.cols(); ++j)
		for (size_t i = 0; i < U.rows(); ++i)
			for (int k = i - 1; k >= 0; --k)
				U(i, j) -= U(k, j) * L(i, k);

	size_t iters = 0, begin = A.cols() * external_i + external_i;
	for (auto i = L._data.begin(); i != L._data.end(); i += L.cols()) {
		std::copy(i, i + L.cols(), A._data.begin() + begin + A.cols() * iters);
		++iters;
	}
	iters = 0;
	for (auto i = U._data.begin(); i != U._data.end(); i += U.cols()) {
		std::copy(i, i + U.cols(), A._data.begin() + begin + A.cols() * iters + external_i);
		++iters;
	}
}

// Substract product of source1 and source2 spans
// from dest span
template <typename T>
void span_substract_product(Matrix<T>& A, const Matrix<T>& L, Matrix<T>& U, size_t external_i) {

	//static_assert(L.rows() != U.rows(), "span_substract_product(): incompatible matrices encountered.");

	for (size_t i = 0; i < A.rows(); ++i)
		for (size_t k = 0; k < L.cols(); ++k)
			for (size_t j = 0; j < A.cols(); ++j)
				A(external_i + i, external_i + j) -= L(external_i + i, k) * U(k, j);
}


//template<typename T>
//void span_set_product_inverseL_by_self(
//	const Matrix<T>& src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
//	Matrix<T>& src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2)
//{
//	if (rows1 != rows2) {
//		std::cout << "inverse_product(): incompatible matrices encountered.";
//		exit(1);
//	}
//
//	for (size_t j = 0; j < cols2; ++j)
//		for (size_t i = 0; i < rows2; ++i)
//			for (int k = i - 1; k >= 0; --k)
//				src2(src2_i + i, src2_j + j) -= src2(src2_i + k, src2_j + j) * src1(src1_i + i, src1_j + k);
//}
//
//// Substract product of source1 and source2 spans
//// from dest span
//template <typename T>
//void span_substract_product(
//	const Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
//	const Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
//	Matrix<T> &dst, size_t dst_i, size_t dst_j) {
//
//	if (cols1 != rows2) throw std::runtime_error("span_substract_product(): incompatible matrices encountered.");
//
//	for (size_t i = 0; i < rows1; ++i)
//		for (size_t k = 0; k < cols1; ++k)
//			for (size_t j = 0; j < cols2; ++j)
//				dst(dst_i + i, dst_j + j) -= src1(src1_i + i, src1_j + k) * src2(src2_i + k, src2_j + j);
//}