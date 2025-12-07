#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include "matrix/sparsematrix.h" 
#include "verfahren/jacobi.h"

using T = double;

// Hilfsfunktion zur Ausgabe einer Matrix (dichtes Format)
void printDenseMatrix(const SparseMatrix<T>& A) {
    std::cout << "\n--- Ausgangsmatrix A (dicht dargestellt) ---\n";
    size_t rows = A.getRows();
    size_t cols = A.getCols();

    std::cout << std::fixed << std::setprecision(4);
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << std::setw(8) << A(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "-------------------------------------------\n";
}

int main() {
    try {
        // --- Setup der Matrix A ---
        const size_t N_ROWS = 3;
        const size_t N_COLS = 3;
        const size_t MAX_ENTRIES_N = 3; 
        
        SparseMatrix<T> A(N_ROWS, N_COLS, MAX_ENTRIES_N);

        // Einträge: [ 10 -1  0 ], [ -1 10 -2 ], [  0 -2 10 ]
        A.addEntry(0, 0, 10.0);
        A.addEntry(0, 1, -1.0);
        A.addEntry(1, 0, -1.0);
        A.addEntry(1, 1, 10.0);
        A.addEntry(1, 2, -2.0);
        A.addEntry(2, 1, -2.0);
        A.addEntry(2, 2, 10.0);
        
        std::vector<T> b = {9.0, 7.0, 8.0};
        std::vector<T> x_start(N_ROWS, 0.0);
        
        const int MAX_ITER = 20;
        const T TOLERANCE = 1e-6;

        // 1. AUSGABE DER MATRIX (wie gewünscht)
        printDenseMatrix(A);
        
        std::cout << "Lösung b: [9.0, 7.0, 8.0]\n";
        std::cout << "Startvektor x^0: [0.0, 0.0, 0.0]\n\n";

        // 2. Aufruf des Jacobi-Lösers
        std::cout << "--- Starte Jacobi-Verfahren ---\n";
        std::vector<T> x_solution = solveJacobi<T>(A, b, x_start, MAX_ITER, TOLERANCE);
        std::cout << "-------------------------------------------\n";

        // 3. Finale Ausgabe der Lösung
        std::cout << "\n** FINALE LÖSUNG **\n";
        std::cout << "Ergebnis x: [";
        for (size_t i = 0; i < N_ROWS; ++i) {
            std::cout << x_solution[i] << (i == N_ROWS - 1 ? "" : ", ");
        }
        std::cout << "]\n";

    } catch (const std::exception& e) {
        std::cerr << "\nFEHLER: " << e.what() << "\n";
    }
    return 0;
}