#ifndef GAUSSSEIDEL_H
#define GAUSSSEIDEL_H

#include "../matrix/sparsematrix.h" 
#include <vector>
#include <stdexcept>

template <typename T>
std::vector<T> solveGaussSeidel(
    const SparseMatrix<T>& A,
    const std::vector<T>& b,
    std::vector<T> x_start,
    int max_iterations,
    T tolerance
);

// Explizite Deklaration der Instanziierung für double
extern template std::vector<double> solveGaussSeidel<double>(
    const SparseMatrix<double>& A,
    const std::vector<double>& b,
    std::vector<double> x_start,
    int max_iterations,
    double tolerance
);

#endif // GAUSSSEIDEL_H