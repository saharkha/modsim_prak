#ifndef JACOBI_H
#define JACOBI_H

#include "../matrix/sparsematrix.h" // Für die Verwendung von SparseMatrix<T>
#include <vector>
#include <stdexcept>

/**
 * @brief Löst das lineare Gleichungssystem A * x = b mit dem Jacobi-Verfahren.
 * * @param A Die Koeffizientenmatrix (SparseMatrix<T>).
 * @param b Der Vektor der rechten Seite.
 * @param x_start Der initiale Lösungsvektor (Startschätzung).
 * @param max_iterations Die maximale Anzahl der Iterationen.
 * @param tolerance Die Abbruchschwelle für die Konvergenzprüfung.
 * @return Der gelöste Vektor x.
 */
template <typename T>
std::vector<T> solveJacobi(
    const SparseMatrix<T>& A,
    const std::vector<T>& b,
    std::vector<T> x_start,
    int max_iterations,
    T tolerance
);

// --- Explizite Deklaration der Instanziierung ---
// Damit der Compiler weiß, dass die Definition für double in der Jacobi.cpp liegt.
extern template std::vector<double> solveJacobi<double>(
    const SparseMatrix<double>& A,
    const std::vector<double>& b,
    std::vector<double> x_start,
    int max_iterations,
    double tolerance
);

#endif // JACOBI_H