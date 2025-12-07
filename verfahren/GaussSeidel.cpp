#include "GaussSeidel.h"
#include <cmath>
#include <iostream>
#include <numeric>

// --- Definition der Template-Funktion ---
template <typename T>
std::vector<T> solveGaussSeidel(
    const SparseMatrix<T>& A,
    const std::vector<T>& b,
    std::vector<T> x_k, // x_k wird in-place überschrieben
    int max_iterations,
    T tolerance
) {
    const size_t n = A.getRows();
    if (b.size() != n || x_k.size() != n) {
        throw std::invalid_argument("Dimension mismatch between matrix and vectors.");
    }
    
    // x_k_old speichert den Zustand VOR der Berechnung der aktuellen Iteration,
    // um den Fehler zu berechnen.
    std::vector<T> x_k_old(n);
    T error = 0.0; 
    
    // 1. Vorbereiten der Diagonalinverse (D^-1)
    std::vector<T> diag_inverse(n);
    for (size_t i = 0; i < n; ++i) {
        T A_ii = A(i, i); 
        if (std::abs(A_ii) < 1e-12) {
            throw std::runtime_error("Diagonal element A_ii is zero (or near zero). Gauss-Seidel method cannot proceed.");
        }
        diag_inverse[i] = 1.0 / A_ii;
    }

    // 2. Haupt-Iterationsschleife
    for (int k = 0; k < max_iterations; ++k) {
        
        // **WICHTIG:** Speichern des Vektors VOR der In-Place-Berechnung
        x_k_old = x_k; 
        
        // Iteration über alle Zeilen (i)
        for (size_t i = 0; i < n; ++i) {
            T sum_terms = 0.0; // Speichert die Summe aller (bekannten) A_ij * x_j Terme
            
            // Verwende den Iterator für die effiziente Summation über die Zeile i
            for (auto it = A.row_begin(i); it != A.row_end(i); ++it) {
                const auto& entry = *it;
                
                if (entry.col_index == A.FILL_INDEX) {
                    break;
                }
                
                size_t j = entry.col_index;
                
                // Nur Off-Diagonalelemente (j != i) tragen zur Summe bei
                if (i != j) {
                    // **DER GAUSS-SEIDEL-UNTERSCHIED:**
                    // Wir verwenden IMMER den aktuellen Zustand von x_k.
                    // Wenn j < i: x_k[j] ist bereits x_j^(k+1).
                    // Wenn j > i: x_k[j] ist immer noch x_j^(k).
                    sum_terms += entry.value * x_k[j]; 
                }
            }
            
            // Gauss-Seidel-Formel: x_i^(k+1) = (1/A_ii) * (b_i - Summe(A_ij * x_j))
            // Die Berechnung wird direkt in x_k[i] gespeichert (in-place)
            x_k[i] = diag_inverse[i] * (b[i] - sum_terms);
        }
        
        // 3. Konvergenzprüfung (L2-Norm)
        T diff_norm_sq = 0.0;
        for (size_t i = 0; i < n; ++i) {
            // Vergleich mit dem Zustand VOR der Iteration
            T diff = x_k[i] - x_k_old[i];
            diff_norm_sq += diff * diff;
        }
        error = std::sqrt(diff_norm_sq); 
        
        if (error < tolerance) {
            std::cout << "Gauss-Seidel konvergiert nach " << k + 1 << " Iterationen.\n";
            return x_k;
        }
        
        // KEIN x_k = x_k_plus_1 Tausch notwendig, da in-place gearbeitet wird!
    }

    std::cout << "Gauss-Seidel erreicht max. Iterationen (" << max_iterations << "). Letzter Fehler: " << error << "\n";
    return x_k;
}

// --- EXPLIZITE INSTANZIIERUNG FÜR DEN TYP DOUBLE ---
template std::vector<double> solveGaussSeidel<double>(
    const SparseMatrix<double>& A,
    const std::vector<double>& b,
    std::vector<double> x_start,
    int max_iterations,
    double tolerance
);