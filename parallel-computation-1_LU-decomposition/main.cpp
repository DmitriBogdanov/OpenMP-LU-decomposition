#include <iostream>
#include <string>

#include "LU_seq.hpp"
#include "static_timer.hpp"


int main(int argc, char *argv[]) {
	const size_t ROWS = 1024;
	const size_t COLS = 1024;
	const size_t BLOCK_SIZE = 8;

	const size_t MIN_SIZE = std::min(ROWS, COLS);
	const size_t SKIP_VERIFICATION_AFTER_SIZE = 1100;

	// Create random matrix
	const auto INITIAL_MATRIX = DMatrix(ROWS, COLS).randomize();

	std::cout << "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n\n";
	std::cout << "INITIAL_MATRIX = " << INITIAL_MATRIX << "\n\n";

	double time1 = 0;
	double time2 = 0;

	// 1) Regular LU decomposition
	{
		auto A = INITIAL_MATRIX;

		// Method
		std::cout << ">>> Regular LU decomposition\n";
		StaticTimer::start();

		LU_seq(A.data(), A.rows(), A.cols());

		time1 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time1 << " sec\n\n";

		// Display
		if (MIN_SIZE < SKIP_VERIFICATION_AFTER_SIZE)
			std::cout
			<< "A = " << A << "\n\n"
			<< "max_element(L * U - INITIAL_MATRIX) = " << verify_LU(A, INITIAL_MATRIX).max_elem() << "\n\n";
		else
			std::cout
			<< "max_element(L * U - INITIAL_MATRIX) = <matrix size exceeds the treshold of " << SKIP_VERIFICATION_AFTER_SIZE << ">\n\n";
	}

	// 2) Block LU decomposition
	if (ROWS == COLS) {
		auto A = INITIAL_MATRIX;

		// Method
		std::cout << ">>> Block LU decomposition\n";
		StaticTimer::start();

		LU_seq_block(A.data(), A.rows(), BLOCK_SIZE);

		time2 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time2 << " sec\n\n";

		// Display
		if (MIN_SIZE < SKIP_VERIFICATION_AFTER_SIZE)
			std::cout
			<< "A = " << A << "\n\n"
			<< "max_element(L * U - INITIAL_MATRIX) = " << verify_LU(A, INITIAL_MATRIX).max_elem() << "\n\n";
		else
			std::cout
			<< "max_element(L * U - INITIAL_MATRIX) = <matrix size exceeds the treshold of " << SKIP_VERIFICATION_AFTER_SIZE << ">\n\n";
	}
	else {
		std::cout << ">>> ROWS != COLS, Block LU decomposition supressed.\n\n";
	}

	std::cout << "Time block / Time seq = " << time2 / time1 << "\n";

	return 0;
}