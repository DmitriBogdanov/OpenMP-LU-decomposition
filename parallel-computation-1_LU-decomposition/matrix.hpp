#pragma once

#include <vector>
#include <algorithm>
#include <ostream>
#include <iomanip>



template <typename T>
struct Matrix {
	// Zero-initialize
	Matrix(size_t rows, size_t cols) : _rows(rows), _cols(cols), _data(rows * cols, 0) {}
	Matrix(size_t size) : _rows(size), _cols(size), _data(size * size, 0) {}

	// Inititalize from data
	Matrix(size_t rows, size_t cols, T *data) : _rows(rows), _cols(cols), _data(data, data + rows * cols) {}
	Matrix(size_t size, T *data) : _rows(size), _cols(size), _data(data, data + size * size) {}
		// Note the use of two iterator constructor 'vector(ptr, ptr + len)'

	// 2D indexation
	T& operator()(size_t i, size_t j) { return _data[i * _cols + j]; }
	const T& operator()(size_t i, size_t j) const { return _data[i * _cols + j]; }

	// 1D indexation
	T& operator[](size_t index) { return _data[index]; }
	const T& operator[](size_t index) const { return _data[index]; }

	// Getters
	size_t rows() const { return _rows; }
	size_t cols() const { return _cols; }

	// Utils
	Matrix<T>& randomize(T min = 0, T max = 100) {
		for (auto &elem : _data) elem = static_cast<T>(min + (max - min) * rand() / (RAND_MAX + 1.));
			// generate random double in [min, max] range and cast to 'T'

		for (size_t k = 0; k < std::min(_rows, _cols); ++k) this->operator()(k, k) *= 5;
			// loosely ensure maximum possible rank

		return *this;
	}

	std::vector<T> _data;
		// direct access for your dirty little needs, use with caution

private:
	size_t _rows;
	size_t _cols;

	
};

using LMatrix = Matrix<long double>;
using DMatrix = Matrix<double>;
using FMatrix = Matrix<float>;
using IMatrix = Matrix<int>;


template<typename T>
std::ostream& operator<< (std::ostream &stream, const Matrix<T> &matrix) {
	for (size_t i = 0; i < matrix.rows(); ++i) {
		for (size_t j = 0; j < matrix.cols(); ++j) stream << std::setw(9) << matrix(i, j);
		stream << '\n';
	}

	return stream;
}