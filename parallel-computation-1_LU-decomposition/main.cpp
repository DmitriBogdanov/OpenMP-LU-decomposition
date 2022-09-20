#include <iostream>

#include "LU_seq.hpp"

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

	/// TEMP (TESTING)
	constexpr size_t ROWS = 8, COLS = 8;

	// Create random matrix
	auto A = DMatrix(ROWS, COLS);
	A.randomize();

	std::cout << "A = \n" << A << "\n\n";

	// LU decomposition
	auto L = DMatrix(ROWS, COLS), U = DMatrix(COLS, COLS);

	//StaticTimer::start();
	//LU_seq(A, L, U);
	//std::cout << "Standart LU decomposition time: " << StaticTimer::end() << "(sec) \n";

	//std::cout << "L = \n" << L << "\n\n" << "U = \n" << U << "\n";
	//std::cout << "L * U = \n" << L * U << "\n";
	
	LU_seq_2(A, L, U);

	/*std::cout << "L = \n" << L << "\n\n" << "U = \n" << U << "\n";
	std::cout << "L * U = \n" << L * U << "\n";*/





	return 0;
}