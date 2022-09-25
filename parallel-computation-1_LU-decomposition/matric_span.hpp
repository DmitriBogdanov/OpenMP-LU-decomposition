#pragma once

#include "matrix.hpp"


// Copy a block into another
template<typename T>
inline void span_copy(
	Matrix<T> &src, size_t src_i, size_t src_j, size_t rows, size_t cols,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

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


template<typename T>
void matrix_mult_wiki_block(const T *A, const T *B, T *C, const int N, const int M, const int K) {
	const int block_size = 8; // 64 = common cache line size
	/*for (int i = 0; i < N; i++) {
		for (int j = 0; j < K; j++) {
			C[K * i + j] = 0;
		}
	}*/
	for (int i0 = 0; i0 < N; i0 += block_size) {
		int imax = i0 + block_size > N ? N : i0 + block_size;

		for (int j0 = 0; j0 < M; j0 += block_size) {
			int jmax = j0 + block_size > M ? M : j0 + block_size;

			for (int k0 = 0; k0 < K; k0 += block_size) {
				int kmax = k0 + block_size > K ? K : k0 + block_size;

				for (int j1 = j0; j1 < jmax; ++j1) {
					int sj = M * j1;

					for (int i1 = i0; i1 < imax; ++i1) {
						int mi = M * i1;
						int ki = K * i1;
						int kij = ki + j1;

						for (int k1 = k0; k1 < kmax; ++k1) {
							C[kij] -= A[mi + k1] * B[sj + k1];
						}
					}
				}
			}
		}
	}
}


template <typename T>
inline void span_substract_product_block(
	const Matrix<T> &src1, size_t src1_i, size_t src1_j, size_t rows1, size_t cols1,
	const Matrix<T> &src2, size_t src2_i, size_t src2_j, size_t rows2, size_t cols2,
	Matrix<T> &dst, size_t dst_i, size_t dst_j) {

	assert(cols1 == rows2 && "span_substract_product(): incompatible matrices encountered.");

	matrix_mult_wiki_block(&(src1._data[0]), &(src2._data[0]), &(dst._data[0]), rows1, cols2, cols1);
}



//void mulBlockHand(double* A, double* B, double* C)
//{
//	double t1, t2;
//	const int bs = 4;
//
//	double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
//	double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
//	double c00, c01, c02, c03, c10, c11, c12, c13, c20, c21, c22, c23, c30, c31, c32, c33;
//
//
//	std::cout << "Proceeding..." << std::endl;
//
//	t1 = clock();
//	for (int bi = 0; bi < n; bi += bs)
//		for (int bj = 0; bj < n; bj += bs)
//		{
//			c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
//			c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
//			c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
//			c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;
//
//			for (int bk = 0; bk < n; bk += bs)
//			{
//				a00 = A[(bi + 0) * n + (bk + 0)]; a01 = A[(bi + 0) * n + (bk + 1)]; a02 = A[(bi + 0) * n + (bk + 2)]; a03 = A[(bi + 0) * n + (bk + 3)];
//				a10 = A[(bi + 1) * n + (bk + 0)]; a11 = A[(bi + 1) * n + (bk + 1)]; a12 = A[(bi + 1) * n + (bk + 2)]; a13 = A[(bi + 1) * n + (bk + 3)];
//				a20 = A[(bi + 2) * n + (bk + 0)]; a21 = A[(bi + 2) * n + (bk + 1)]; a22 = A[(bi + 2) * n + (bk + 2)]; a23 = A[(bi + 2) * n + (bk + 3)];
//				a30 = A[(bi + 3) * n + (bk + 0)]; a31 = A[(bi + 3) * n + (bk + 1)]; a32 = A[(bi + 3) * n + (bk + 2)]; a33 = A[(bi + 3) * n + (bk + 3)];
//
//				b00 = B[(bj + 0) * n + (bk + 0)]; b01 = B[(bj + 0) * n + (bk + 1)]; b02 = B[(bj + 0) * n + (bk + 2)]; b03 = B[(bj + 0) * n + (bk + 3)];
//				b10 = B[(bj + 1) * n + (bk + 0)]; b11 = B[(bj + 1) * n + (bk + 1)]; b12 = B[(bj + 1) * n + (bk + 2)]; b13 = B[(bj + 1) * n + (bk + 3)];
//				b20 = B[(bj + 2) * n + (bk + 0)]; b21 = B[(bj + 2) * n + (bk + 1)]; b22 = B[(bj + 2) * n + (bk + 2)]; b23 = B[(bj + 2) * n + (bk + 3)];
//				b30 = B[(bj + 3) * n + (bk + 0)]; b31 = B[(bj + 3) * n + (bk + 1)]; b32 = B[(bj + 3) * n + (bk + 2)]; b33 = B[(bj + 3) * n + (bk + 3)];
//
//				c00 += a00 * b00 + a01 * b01 + a02 * b02 + a03 * b03;
//				c01 += a00 * b10 + a01 * b11 + a02 * b12 + a03 * b13;
//				c02 += a00 * b20 + a01 * b21 + a02 * b22 + a03 * b23;
//				c03 += a00 * b30 + a01 * b31 + a02 * b32 + a03 * b33;
//
//				c10 += a10 * b00 + a11 * b01 + a12 * b02 + a13 * b03;
//				c11 += a10 * b10 + a11 * b11 + a12 * b12 + a13 * b13;
//				c12 += a10 * b20 + a11 * b21 + a12 * b22 + a13 * b23;
//				c13 += a10 * b30 + a11 * b31 + a12 * b32 + a13 * b33;
//
//				c20 += a20 * b00 + a21 * b01 + a22 * b02 + a23 * b03;
//				c21 += a20 * b10 + a21 * b11 + a22 * b12 + a23 * b13;
//				c22 += a20 * b20 + a21 * b21 + a22 * b22 + a23 * b23;
//				c23 += a20 * b30 + a21 * b31 + a22 * b32 + a23 * b33;
//
//				c30 += a30 * b00 + a31 * b01 + a32 * b02 + a33 * b03;
//				c31 += a30 * b10 + a31 * b11 + a32 * b12 + a33 * b13;
//				c32 += a30 * b20 + a31 * b21 + a32 * b22 + a33 * b23;
//				c33 += a30 * b30 + a31 * b31 + a32 * b32 + a33 * b33;
//			}
//
//			C[(bi + 0) * n + (bj + 0)] = c00; C[(bi + 0) * n + (bj + 1)] = c01; C[(bi + 0) * n + (bj + 2)] = c02; C[(bi + 0) * n + (bj + 3)] = c03;
//			C[(bi + 1) * n + (bj + 0)] = c10; C[(bi + 1) * n + (bj + 1)] = c11; C[(bi + 1) * n + (bj + 2)] = c12; C[(bi + 1) * n + (bj + 3)] = c13;
//			C[(bi + 2) * n + (bj + 0)] = c20; C[(bi + 2) * n + (bj + 1)] = c21; C[(bi + 2) * n + (bj + 2)] = c22; C[(bi + 2) * n + (bj + 3)] = c23;
//			C[(bi + 3) * n + (bj + 0)] = c30; C[(bi + 3) * n + (bj + 1)] = c31; C[(bi + 3) * n + (bj + 2)] = c32; C[(bi + 3) * n + (bj + 3)] = c33;
//		}
//	t2 = clock();
//
//	std::cout << "Time BlockHand = " << (double)(t2 - t1) / CLOCKS_PER_SEC << " sec" << std::endl;
//	std::cout << "C[n/2][n/2] = " << C[n / 2 * n + n / 2] << std::endl;
//
//}