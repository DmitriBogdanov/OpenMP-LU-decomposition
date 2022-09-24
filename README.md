# Parallel Computation / LU Decomposition

Contains sequential and parallelized implementations of following algorithms:

* LU decomposition
* Block LU decomposition

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: -TODO-
* Requires C++17 support

## Usage

-TODO-

## Version history
* 00.07
    * Fixed incorrect LU split for rectangular matrices with more rows

* 00.06
    * Implemented BLAS3 LU factorization

* 00.05
    * Bugfixes in BLAS3 decomposition, narrowed down issues to step (3)

* 00.04
    * Refactored BLAS2 LU decomposition
    * Implemented operations for matrix spans
    * Partially implemented LU decomposition BLAS3 for blocks of set size (b = 2)
    * General refactors

* 00.03
    * Implemented sequential LU decomposition BLAS2 in place

* 00.02
    * Created suitable matrix implementation
    * Implemented sequential LU decomposition BLAS2

* 00.01
    * Implemented thread safe timing

## License

This project is licensed under the MIT License - see the LICENSE.md file for details