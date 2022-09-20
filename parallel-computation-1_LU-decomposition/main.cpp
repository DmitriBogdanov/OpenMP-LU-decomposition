#include <iostream>

#include "static_timer.hpp"
#include "LU_seq.h"



int main(int argc, char *argv[]) {
	/// TEMP (TESTING)
	constexpr size_t ROWS = 6, COLS = 6;

	// Create random matrix
	auto A = DMatrix(ROWS, COLS);
	A.randomize();

	std::cout << A << "\n\n";

	// LU decomposition
	auto L = DMatrix(ROWS, COLS), U = DMatrix(ROWS, COLS);

	LU_seq(A, L, U);
	

	std::cout << L << "\n\n" << U << "\n";

	return 0;
}