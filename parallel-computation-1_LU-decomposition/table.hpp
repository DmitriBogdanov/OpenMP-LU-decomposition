#pragma once

#include <iostream>
#include <iomanip>
#include <string>


template<typename Object>
void table_add_1(const Object &obj) {
	std::cout
		<< " | "
		<< std::setw(17) << obj
		<< " | ";
}


template<typename Object>
void table_add_2(const Object &obj) {
	std::cout
		<< std::setw(12) << obj
		<< " | ";
}


template<typename Object>
void table_add_3(const Object &obj) {
	std::cout
		<< std::setw(12) << obj
		<< " | ";
}


template<typename Object>
void table_add_4(const Object &obj) {
	std::cout
		<< std::setw(12) << obj
		<< " |\n";
}


void table_hline() {
	std::string str1("");
	std::string str2("");
	std::string str3("");
	std::string str4("");
	str1.insert(0, 1 + 17 + 1, '-');
	str2.insert(0, 1 + 12 + 1, '-');
	str3.insert(0, 1 + 12 + 1, '-');
	str4.insert(0, 1 + 12 + 1, '-');

	std::cout << " |" << str1 << "|" << str2 << "|" << str3 << "|" << str4 << "|\n";
}