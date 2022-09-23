#include <iostream>

#include "LU_seq.hpp"
#include "static_timer.hpp"

void cpp_standart()
{
	if (__cplusplus == 201703L) std::cout << "C++17\n";
	else if (__cplusplus == 201402L) std::cout << "C++14\n";
	else if (__cplusplus == 201103L) std::cout << "C++11\n";
	else if (__cplusplus == 199711L) std::cout << "C++98\n";
	else std::cout << "pre-standard C++\n";
}


int main(int argc, char *argv[]) {
	cpp_standart();

	constexpr size_t ROWS = 3000;
	constexpr size_t COLS = 3000;

	// Create random matrix
	const auto INITIAL_MATRIX = DMatrix(ROWS, COLS).randomize();

	std::cout << "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n";
	std::cout << "INITIAL_MATRIX =\n" << INITIAL_MATRIX << "\n\n";
	double time1, time2;
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
		std::cout << "Completed in " << time1 << "sec\n";

		// Display
		split_into_LU(A, L, U);
		std::cout
			<< "A = \n" << A << "\n\n"
			<< "L = \n" << L << "\n\n"
			<< "U = \n" << U << "\n\n"
			<< "max_norm(L * U - INITIAL_MATRIX) = \n" << (L * U - INITIAL_MATRIX).max_elem() << "\n";
	}

	// 2) Block LU decomposition
	{
		auto A = INITIAL_MATRIX;
		auto L = DMatrix(ROWS, std::min(ROWS, COLS));
		auto U = DMatrix(std::min(ROWS, COLS), COLS);

		// Method
		std::cout << ">>> Block LU decomposition\n";
		StaticTimer::start();
		const size_t b = 4;
		LU_seq_block(A, b);

		time2 = StaticTimer::end() / 1000.;
		std::cout << "Completed in " << time2 << "sec\n";

		// Display
		split_into_LU(A, L, U);
		std::cout
			<< "A = \n" << A << "\n\n"
			<< "L = \n" << L << "\n\n"
			<< "U = \n" << U << "\n\n"
			<< "max_norm(L * U - INITIAL_MATRIX) = \n" << (L * U - INITIAL_MATRIX).max_elem() << "\n";
	}
	std::cout << "Time block / Time seq = " << time2 / time1 << "\n";

	return 0;
}