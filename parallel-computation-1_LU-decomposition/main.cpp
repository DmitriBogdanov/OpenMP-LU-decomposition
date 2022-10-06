#include <iostream>

#include "LU_seq.hpp"
#include "static_timer.hpp"


int main(int argc, char *argv[]) {
	const size_t ROWS = pow(2, 11);
	const size_t COLS = pow(2, 11);
	const size_t BLOCK_SIZE = pow(2, 5);

	const size_t num_of_experiments = 1;
	std::vector<double> ratios(num_of_experiments, 0.);
	DMatrix A(ROWS, COLS);
	DMatrix INITIAL_MATRIX(ROWS, COLS);

	for (size_t i = 0; i < num_of_experiments; ++i)
	{
		std::cout << "EXPERIMENT ¹" << i + 1 << "--------------------START--------------------\n";
		// Create random matrix
		INITIAL_MATRIX.randomize();

		std::cout << "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n\n";
		std::cout << "BLOCK_DIMENSIONS = (" << BLOCK_SIZE << ", " << BLOCK_SIZE << ")\n\n";
		std::cout << "INITIAL_MATRIX = " << INITIAL_MATRIX << "\n\n";

		double time1 = 0, time2 = 0;

		// 1) Regular LU decomposition
		{
			A = INITIAL_MATRIX;
			/*	auto L = DMatrix(ROWS, std::min(ROWS, COLS));
				auto U = DMatrix(std::min(ROWS, COLS), COLS);*/

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
			A = INITIAL_MATRIX;
			/*auto L = DMatrix(ROWS, std::min(ROWS, COLS));
			auto U = DMatrix(std::min(ROWS, COLS), COLS);*/

			// Method
			std::cout << ">>> Block LU decomposition\n";
			StaticTimer::start();

			LU_seq_block(A, BLOCK_SIZE);

			time2 = StaticTimer::end() / 1000.;
			std::cout << "Completed in " << time2 << "sec\n\n";

			// Display
		/*	split_into_LU(A, L, U);*/
		/*	std::cout
				<< "A = " << A << "\n\n"
				<< "max_norm(L * U - INITIAL_MATRIX) = " << (product_LU(A) - INITIAL_MATRIX).max_elem() << "\n\n";*/
		}
		else {
			std::cout << ">>> ROWS != COLS, Block LU decomposition supressed.\n\n";
		}
		ratios[i] = time2 / time1;
		std::cout << "Time block / Time seq = " << time2 / time1 << "\n";
		std::cout << "EXPERIMENT ¹" << i + 1 << "--------------------FINISH--------------------\n";
	}
	double mean = 0;
	for (size_t i = 0; i < ratios.size(); ++i)
	{
		mean += ratios[i];
	}
	std::cout << "Mean value of ( num of experiments = " << ratios.size() << " ) Time block / Time seq = " << mean / ratios.size() << "\n";


	return 0;
}