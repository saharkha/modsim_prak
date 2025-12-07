#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm> // Für std::copy

// WICHTIG: Korrekte Includes für Ihre Ordnerstruktur
#include "matrix/sparsematrix.h" 
#include "verfahren/jacobi.h"
#include "verfahren/GaussSeidel.h" // NEU: Inclusion des Gauss-Seidel Lösers

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

// Hilfsfunktion zur Ausgabe des Ergebnisvektors
void printSolution(const std::vector<T>& x_solution) {
    size_t n = x_solution.size();
    std::cout << "Ergebnis x: [";
    for (size_t i = 0; i < n; ++i) {
        std::cout << std::fixed << std::setprecision(6) << x_solution[i] << (i == n - 1 ? "" : ", ");
    }
    std::cout << "]\n";
}

int main() {
    try {
        // --- Setup der Matrix A (Diagonaldmoninant für Konvergenz) ---
        const size_t N_ROWS = 3;
        const size_t N_COLS = 3;
        const size_t MAX_ENTRIES_N = 3; 
        
        SparseMatrix<T> A(N_ROWS, N_COLS, MAX_ENTRIES_N);

        // Matrix A: [ 10 -1  0 ], [ -1 10 -2 ], [  0 -2 10 ]
        A.addEntry(0, 0, 10.0);
        A.addEntry(0, 1, -1.0);
        A.addEntry(1, 0, -1.0);
        A.addEntry(1, 1, 10.0);
        A.addEntry(1, 2, -2.0);
        A.addEntry(2, 1, -2.0);
        A.addEntry(2, 2, 10.0);
        
        // Vektor b: Für die bekannte Lösung x = [1, 1, 1]
        std::vector<T> b = {9.0, 7.0, 8.0};
        
        // Parameter
        const int MAX_ITER = 50;
        const T TOLERANCE = 1e-6;

        // 1. Ausgabe der Matrix
        printDenseMatrix(A);
        
        std::cout << "Lösung b: [9.0, 7.0, 8.0]\n";
        std::cout << "Gesuchte Lösung: x = [1.0, 1.0, 1.0]\n\n";

        // --- 2. Jacobi-Verfahren ---
        std::cout << "===========================================\n";
        std::cout << "          STARTE JACOBI-VERFAHREN          \n";
        std::cout << "===========================================\n";
        
        // Wichtig: x_start wird an solveJacobi als Kopie übergeben
        std::vector<T> x_start(N_ROWS, 0.0); 
        std::vector<T> x_jacobi = solveJacobi<T>(A, b, x_start, MAX_ITER, TOLERANCE);
        
        std::cout << "Jacobi FINALE LÖSUNG:\n";
        printSolution(x_jacobi);
        std::cout << "-------------------------------------------\n";

        // --- 3. Gauss-Seidel-Verfahren ---
        std::cout << "===========================================\n";
        std::cout << "        STARTE GAUSS-SEIDEL-VERFAHREN      \n";
        std::cout << "===========================================\n";

        // Wichtig: Erneute Startschätzung für den fairen Vergleich
        std::vector<T> x_start_gs(N_ROWS, 0.0);
        std::vector<T> x_gauss_seidel = solveGaussSeidel<T>(A, b, x_start_gs, MAX_ITER, TOLERANCE);

        std::cout << "Gauss-Seidel FINALE LÖSUNG:\n";
        printSolution(x_gauss_seidel);
        std::cout << "-------------------------------------------\n";


    } catch (const std::exception& e) {
        std::cerr << "\nFEHLER: " << e.what() << "\n";
    }
    return 0;
}