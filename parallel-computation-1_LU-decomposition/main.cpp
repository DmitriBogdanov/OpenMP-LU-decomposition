#include <iostream>

#include "LU_seq.hpp"
#include "static_timer.hpp"


int main(int argc, char *argv[]) {
	const size_t ROWS = 256;
	const size_t COLS = 256;
	const size_t BLOCK_SIZE = 16;

	// Create random matrix
	const auto INITIAL_MATRIX = DMatrix(ROWS, COLS).randomize();

	std::cout << "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n\n";
	std::cout << "INITIAL_MATRIX = " << INITIAL_MATRIX << "\n\n";

	double time1 = 0, time2 = 0;

	// 1) Regular LU decomposition
	{
		auto A = INITIAL_MATRIX;
		auto L = DMatrix(ROWS, std::min(ROWS, COLS));
		auto U = DMatrix(std::min(ROWS, COLS), COLS);

		// Method
		std::cout << ">>> Regular LU decomposition\n";
		StaticTimer::start();

		LU_seq(A);

		time1 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time1 << "sec\n\n";

		// Display
		/*split_into_LU(A, L, U);
		std::cout
			<< "A = " << A << "\n\n"
			<< "max_norm(L * U - INITIAL_MATRIX) = " << (L * U - INITIAL_MATRIX).max_elem() << "\n\n";*/
	}

	// 2) Block LU decomposition
	if (ROWS == COLS) {
		auto A = INITIAL_MATRIX;
		auto L = DMatrix(ROWS, std::min(ROWS, COLS));
		auto U = DMatrix(std::min(ROWS, COLS), COLS);

		// Method
		std::cout << ">>> Block LU decomposition\n";
		StaticTimer::start();

		LU_seq_block(A, BLOCK_SIZE);

		time2 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time2 << "sec\n\n";

		// Display
		split_into_LU(A, L, U);
		std::cout
			<< "A = " << A << "\n\n"
			<< "max_norm(L * U - INITIAL_MATRIX) = " << (L * U - INITIAL_MATRIX).max_elem() << "\n\n"
			<< "max_norm(L * U - INITIAL_MATRIX) = " << (product_LU(A) - INITIAL_MATRIX).max_elem() << "\n\n";
	}
	else {
		std::cout << ">>> ROWS != COLS, Block LU decomposition supressed.\n\n";
	}

	std::cout << "Time block / Time seq = " << time2 / time1 << "\n";

	return 0;
}