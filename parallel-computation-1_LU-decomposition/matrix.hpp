#pragma once

#include <vector>
#include <algorithm>
#include <ostream>
#include <iomanip>

template <typename T>
struct Matrix {
	// Zero-initialize
	Matrix(size_t rows, size_t cols, T var = static_cast<T>(0)) : _rows(rows), _cols(cols), _data(rows * cols, var) {}
	Matrix(size_t size, T var = static_cast<T>(0)) : _rows(size), _cols(size), _data(size * size, var) {}

	// Initialize from data
	Matrix(size_t rows, size_t cols, T *data) : _rows(rows), _cols(cols), _data(data, data + rows * cols) {}
	Matrix(size_t size, T *data) : _rows(size), _cols(size), _data(data, data + size * size) {}
		// Note the use of two iterator constructor 'vector(ptr, ptr + len)'

	// Initialize from Matrix itself
	template <typename every_type>
	Matrix(Matrix<every_type> matrix) : _rows(matrix.rows()), _cols(matrix.cols()), _data(*matrix._data) {}

	// 2D indexation
	T& operator()(size_t i, size_t j) { return _data[i * _cols + j]; }
	const T& operator()(size_t i, size_t j) const { return _data[i * _cols + j]; }

	// 1D indexation
	T& operator[](size_t index) { return _data[index]; }
	const T& operator[](size_t index) const { return _data[index]; }

	// Getters
	size_t rows() const { return _rows; }
	size_t cols() const { return _cols; }

	//Standart matrix multiplication
	Matrix<T> operator*(const Matrix<T> &other) {
		if (this->cols() != other.rows()) {
			std::cout << "operator*(): incompatible matrices encountered.\n";
			exit(1);
		};

		Matrix<T> res(this->rows(), other.cols());

		for (size_t i = 0; i < this->rows(); ++i)
			for (size_t j = 0; j < other.cols(); ++j)
				for (size_t k = 0; k < this->cols(); ++k)
					res(i, j) += this->operator()(i, k) * other(k, j);
		return res;
	}

	Matrix<T> operator-(const Matrix<T> &other) {
		if (this->rows() != other.rows() || this->cols() != other.cols()) {
			std::cout << "operator-(): incompatible matrices encountered.\n";
			exit(1);
		};

		Matrix<T> res = *this;

		for (size_t i = 0; i < this->rows(); ++i)
			for (size_t j = 0; j < this->cols(); ++j)
					res(i, j) -= other(i, j);
		return res;
	}

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
std::ostream& operator<<(std::ostream &stream, const Matrix<T> &matrix) {
	constexpr size_t MAX_DISPLAYED_ROWS = 12;
	constexpr size_t MAX_DISPLAYED_COLS = 7;

	if (matrix.rows() <= MAX_DISPLAYED_ROWS && matrix.cols() <= MAX_DISPLAYED_COLS) {
		for (size_t i = 0; i < matrix.rows(); ++i) {
			stream << std::setw(4) << '[';
			for (size_t j = 0; j < matrix.cols(); ++j) 
				if (std::abs(matrix(i, j)) < 1000) stream << std::setw(10) << std::setprecision(4) << std::defaultfloat << matrix(i, j);
				else stream << std::setw(10) << std::setprecision(1) << std::scientific << matrix(i, j);
			stream << std::setw(4) << ']' << '\n';
		}
	}
	else {
		stream << "<supressed>";
	}

	return stream;
}