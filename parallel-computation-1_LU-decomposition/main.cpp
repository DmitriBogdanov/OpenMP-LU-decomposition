#include <iostream>

#include "LU_seq.hpp"
#include "static_timer.hpp"

void cpp_standart()
{
	if (_MSVC_LANG == 201703L) std::cout << "C++17\n";
	else if (_MSVC_LANG == 201402L) std::cout << "C++14\n";
	else if (_MSVC_LANG == 201103L) std::cout << "C++11\n";
	else if (_MSVC_LANG == 199711L) std::cout << "C++98\n";
	else std::cout << "pre-standard C++\n\n";
}


int main(int argc, char *argv[]) {
	cpp_standart();

	const size_t ROWS = 10;
	const size_t COLS = 10;
	const size_t BLOCK_SIZE = 2;

	// Create random matrix
	const auto INITIAL_MATRIX = DMatrix(ROWS, COLS).randomize();

	std::cout << "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n\n";
	std::cout << "INITIAL_MATRIX = " << INITIAL_MATRIX << "\n\n";

	double time1 = 0, time2 = 0;

	// 1) Regular LU decomposition
	{
		DMatrix A(INITIAL_MATRIX);

		// Method
		std::cout << ">>> Regular LU decomposition\n";
		StaticTimer::start();

		LU_seq(A);

		time1 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time1 << "sec\n\n";

		std::cout
			<< "A = " << A << "\n\n"
			<< "L * U = " << product_LU(A) << "\n\n";
	}

	// 2) Block LU decomposition
	if (ROWS == COLS) {
		DMatrix A(INITIAL_MATRIX);

		// Method
		std::cout << ">>> Block LU decomposition\n";
		StaticTimer::start();

		LU_seq_block(A, BLOCK_SIZE);

		time2 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time2 << "sec\n\n";

		std::cout
			<< "A = " << A << "\n\n"
			<< "L * U = " << product_LU(A) << "\n\n";
	}
	else {
		std::cout << ">>> ROWS != COLS, Block LU decomposition supressed.\n\n";
	}

	std::cout << "Time block / Time seq = " << time2 / time1 << "\n";

	return 0;
}