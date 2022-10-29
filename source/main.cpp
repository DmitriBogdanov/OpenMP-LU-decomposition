#include "matrix.hpp"
#include "LU_serial.hpp"
#include "LU_parallel.hpp"
#include "static_timer.hpp"
#include "table.hpp"


// Config
const size_t ROWS = 1024 * 2;
const size_t COLS = 1024 * 2;
const size_t BLOCK_SIZE = 32;
const int THREAD_CAP = 4; // to limit threads below 'omp_get_max_threads()'


int main(int argc, char *argv[]) {
	// Consts
	const size_t MAX_SIZE = std::max(ROWS, COLS);
	const size_t SKIP_VERIFICATION_AFTER_SIZE = 1100;
		// since solution verification sometimes takes longer than
		// the method itself, we skip it for large enough sizes


	// Set up OpenMP
	const int MAX_THREADS = omp_get_max_threads();

	omp_set_num_threads(std::min(THREAD_CAP, MAX_THREADS));

	// Create random matrix
	const auto INITIAL_MATRIX = DMatrix(ROWS, COLS).randomize();

	std::cout
		<< "MATRIX_DIMENSIONS = (" << ROWS << ", " << COLS << ")\n"
		<< "BLOCK_SIZE = " << BLOCK_SIZE << "\n\n"
		<< "Using " << std::min(THREAD_CAP, MAX_THREADS) << " threads out of " << MAX_THREADS << "\n\n";

	// System messages
	if (ROWS != COLS)
		std::cout << "> NOTE: ROWS != COLS, Block LU decomposition supressed.\n\n";
	if (ROWS % BLOCK_SIZE != 0)
		std::cout << "> NOTE: ROWS is not a multiple of BLOCK_SIZE, Block LU decomposition supressed.\n\n";
	if (MAX_SIZE > SKIP_VERIFICATION_AFTER_SIZE)
		std::cout << "> NOTE: Matrix size exceeds " << SKIP_VERIFICATION_AFTER_SIZE << ", verification supressed\n\n";

	// Draw a table
	table_add_1("Method");
	table_add_2("Time (sec)");
	table_add_3("Error");
	table_add_4("Speedup");
	table_hline();

	double regularSerialTime = -1;
	double regularParallelTime = -1;
	double blockSerialTime = -1;
	double blockParallelTime = -1;

	// 1) LU
	{
		auto A = INITIAL_MATRIX;

		// 1. Method
		table_add_1("LU");

		// 2. Time
		StaticTimer::start();
		LU_serial(A.data(), A.rows(), A.cols());
		regularSerialTime = StaticTimer::end();

		table_add_2(regularSerialTime);

		// 3. Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// 4. Speedup
		table_add_4(regularSerialTime / regularSerialTime);
	}

	// 2) Parallel LU
	{
		auto A = INITIAL_MATRIX;

		// 1. Method
		table_add_1("Parallel LU");

		// 2. Time
		StaticTimer::start();
		LU_parallel(A.data(), A.rows(), A.cols());
		regularParallelTime = StaticTimer::end();

		table_add_2(regularParallelTime);

		// 3. Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// 4. Speedup
		table_add_4(regularSerialTime / regularParallelTime);
	}

	// 3) Block LU
	if (ROWS == COLS && ROWS % BLOCK_SIZE == 0) {
		auto A = INITIAL_MATRIX;

		// 1. Method
		table_add_1("Block LU");

		// 2. Time
		StaticTimer::start();
		blockLU_serial(A.data(), A.rows(), BLOCK_SIZE);
		blockSerialTime = StaticTimer::end();

		table_add_2(blockSerialTime);

		// 3. Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// 4. Speedup
		table_add_4(regularSerialTime / blockSerialTime);
	}

	// 4) Parallel Block LU
	if (ROWS == COLS && ROWS % BLOCK_SIZE == 0) {
		auto A = INITIAL_MATRIX;

		// 1. Method
		table_add_1("Parallel Block LU");

		// 2. Time
		StaticTimer::start();
		blockLU_parallel(A.data(), A.rows(), BLOCK_SIZE);
		blockParallelTime = StaticTimer::end();

		table_add_2(blockParallelTime);

		// 3. Err
		if (MAX_SIZE < SKIP_VERIFICATION_AFTER_SIZE) table_add_3(verify_LU(A, INITIAL_MATRIX));
		else table_add_3("<supressed>");

		// 4. Speedup
		table_add_4(blockSerialTime / blockParallelTime);
	}

	return 0;
}