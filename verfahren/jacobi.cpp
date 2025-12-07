#include "jacobi.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

// --- Definition der Template-Funktion ---
template <typename T>
std::vector<T> solveJacobi(
    const SparseMatrix<T>& A,
    const std::vector<T>& b,
    std::vector<T> x_start,
    int max_iterations,
    T tolerance
) {
    const size_t n = A.getRows();
    if (b.size() != n || x_start.size() != n) {
        throw std::invalid_argument("Dimension mismatch between matrix and vectors.");
    }
    
    std::vector<T> x_k = x_start;
    std::vector<T> x_k_plus_1(n);
    
    // Wichtig: 'error' VOR der Schleife deklarieren, damit sie am Ende sichtbar ist.
    T error = 0.0; 
    
    // 1. Vorbereiten der Diagonalinverse (D^-1)
    std::vector<T> diag_inverse(n);
    for (size_t i = 0; i < n; ++i) {
        // Direkter Zugriff auf Diagonalelement, langsam, aber für D^-1 notwendig.
        T A_ii = A(i, i); 
        if (std::abs(A_ii) < 1e-12) {
            throw std::runtime_error("Diagonal element A_ii is zero (or near zero). Jacobi method cannot proceed.");
        }
        diag_inverse[i] = 1.0 / A_ii;
    }

    // 2. Haupt-Iterationsschleife
    for (int k = 0; k < max_iterations; ++k) {
        
        // Berechnung des neuen Vektors x^(k+1)
        for (size_t i = 0; i < n; ++i) {
            T off_diag_sum = 0.0; 
            
            // Verwende den Iterator für die effiziente Summation über die Zeile i
            for (auto it = A.row_begin(i); it != A.row_end(i); ++it) {
                const auto& entry = *it;
                
                if (entry.col_index == A.FILL_INDEX) {
                    break;
                }
                
                size_t j = entry.col_index;
                
                // Nur Off-Diagonalelemente (j != i) tragen zur Summe bei
                if (i != j) {
                    off_diag_sum += entry.value * x_k[j]; 
                }
            }
            
            // Jacobi-Formel: x_i^(k+1) = (1/A_ii) * (b_i - Summe(A_ij * x_j^(k)))
            x_k_plus_1[i] = diag_inverse[i] * (b[i] - off_diag_sum);
        }
        
        // 3. Konvergenzprüfung (L2-Norm)
        T diff_norm_sq = 0.0;
        for (size_t i = 0; i < n; ++i) {
            T diff = x_k_plus_1[i] - x_k[i];
            diff_norm_sq += diff * diff;
        }
        error = std::sqrt(diff_norm_sq); 
        
        if (error < tolerance) {
            std::cout << "Jacobi konvergiert nach " << k + 1 << " Iterationen.\n";
            return x_k_plus_1;
        }
        
        // 4. Vorbereitung für die nächste Iteration: x_k = x_k_plus_1
        x_k = x_k_plus_1; 
    }

    // Konsolenausgabe, wenn maximale Iterationen erreicht wurden.
    // 'error' ist hier sichtbar und enthält den letzten berechneten Fehler.
    std::cout << "Jacobi erreicht max. Iterationen (" << max_iterations << "). Letzter Fehler: " << error << "\n";
    return x_k_plus_1;
}

// --- EXPLIZITE INSTANZIIERUNG FÜR DEN TYP DOUBLE ---
template std::vector<double> solveJacobi<double>(
    const SparseMatrix<double>& A,
    const std::vector<double>& b,
    std::vector<double> x_start,
    int max_iterations,
    double tolerance
);