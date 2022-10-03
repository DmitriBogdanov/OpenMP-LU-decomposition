#include <iostream>
#include <string>

#include "matrix.hpp"
#include "LU_seq.hpp"
#include "LU_par.hpp"
#include "static_timer.hpp"
#include "table.hpp"


int main(int argc, char *argv[]) {
	const size_t ROWS = 1024;
	const size_t COLS = 1024;
	const size_t BLOCK_SIZE = 32;

	const size_t MAX_SIZE = std::max(ROWS, COLS);
	const size_t SKIP_VERIFICATION_AFTER_SIZE = 1100;

	// Create random matrix
	const auto INITIAL_MATRIX = DMatrix(ROWS, COLS).randomize();

	std::cout
		<< "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n"
		<< "BLOCK_SIZE = " << BLOCK_SIZE << "\n\n"
		<< "INITIAL_MATRIX = " << INITIAL_MATRIX << "\n\n";

	double defaultTime = 0;
	double currentTime = 0;

	// Display
	if (ROWS != COLS)
		std::cout << "> NOTE: ROWS != COLS, Block LU decomposition supressed.\n";
	if (MAX_SIZE > SKIP_VERIFICATION_AFTER_SIZE)
		std::cout << "> NOTE: Matrix size exceeds " << SKIP_VERIFICATION_AFTER_SIZE << ", verification supressed\n";
	std::cout << "\n";

	table_add_1("Method");
	table_add_2("Time (sec)");
	table_add_3("Error");
	table_add_4("Speedup");
	table_hline();

	// 1) LU
	{
		auto A = INITIAL_MATRIX;

		// Method
		table_add_1("LU");

		// Time
		StaticTimer::start();

		LU_seq(A.data(), A.rows(), A.cols());

		defaultTime = StaticTimer::end() / 1000.;
		currentTime = defaultTime;

		table_add_2(currentTime);

		// Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// Speedup
		table_add_4(currentTime / defaultTime);
	}

	// 2) Block LU
	if (ROWS == COLS) {
		auto A = INITIAL_MATRIX;

		// Method
		table_add_1("Block LU");

		// Time
		StaticTimer::start();

		LU_seq_block(A.data(), A.rows(), BLOCK_SIZE);

		currentTime = StaticTimer::end() / 1000.;

		table_add_2(currentTime);

		// Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// Speedup
		table_add_4(currentTime / defaultTime);
	}

	// 3) Parallel LU
	{
		auto A = INITIAL_MATRIX;

		// Method
		table_add_1("Parallel LU");

		// Time
		StaticTimer::start();

		LU_par(A.data(), A.rows(), A.cols());

		currentTime = StaticTimer::end() / 1000.;

		table_add_2(currentTime);

		// Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// Speedup
		table_add_4(currentTime / defaultTime);
	}

	// 4) Parallel Block LU
	if (ROWS == COLS) {
		auto A = INITIAL_MATRIX;

		// Method
		table_add_1("Parallel Block LU");

		// Time
		StaticTimer::start();

		LU_par_block(A.data(), A.rows(), BLOCK_SIZE);

		currentTime = StaticTimer::end() / 1000.;

		table_add_2(currentTime);

		// Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// Speedup
		table_add_4(currentTime / defaultTime);
	}

	return 0;
}