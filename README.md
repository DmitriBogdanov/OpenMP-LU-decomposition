# OpenMP LU Decomposition

Contains serial and parallel (OpenMP) implementations of following algorithms:

* LU decomposition
* Blocked LU decomposition

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: Intel C++ Compiler
* Requires C++17 support
* Requires OpenMP support

## Usage

Adjust in ROWS, COLS and BLOCK_SIZE in "main.cpp" to configure testing parameters. Block decomposition assumes square matrix with size being a multiple of BLOCK_SIZE.

## Version history

* 01.01
    * Code cleanup and renames across the board
    * Fixed typos in names, comments and docs
    * Added multiple comments with full method descriptions and notes about technical implementation
    * Added safety guards for sizes that are not a multiple of BLOCK_SIZE

* 00.11
    * Removed Windows-specific calls to allow compilation on Linux clusters
    * Optimized  parallel block LU

* 00.10
    * Added date of the computation to output
    * Implemented parallel block LU

* 00.09
    * Changed output style to a table
    * Implemented parallel LU using OpenMP
    * Optimized verify_LU() to avoid allocations

* 00.08
    * Reimplemented LU decomposition and corresponding span methods using C-style arrays and pointer arithmetic
    * Optimized cache use in blocked method
    * Added optimized method for verifying decomposition corectness

* 00.07
    * Fixed incorrect LU split for rectangular matrices with more rows

* 00.06
    * Implemented functional blocked LU factorization

* 00.05
    * Bugfixes in blocked decomposition, narrowed down issues to step (3)

* 00.04
    * Refactored LU decomposition
    * Implemented operations for matrix spans
    * Partially implemented LU decomposition for blocks of set size (b = 2)
    * General refactors

* 00.03
    * Implemented serial LU decomposition in place

* 00.02
    * Created suitable matrix implementation
    * Implemented sequential LU decomposition

* 00.01
    * Implemented thread safe timing

## License

This project is licensed under the MIT License - see the LICENSE.md file for details