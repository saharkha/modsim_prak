#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm> // Für std::copy
#include <cmath>
#include <string>

// WICHTIG: Korrekte Includes für Ordnerstruktur
#include "matrix/sparsematrix.h" 
#include "verfahren/jacobi.h"
#include "verfahren/GaussSeidel.h" 
#include "Mehrgitterverfahren/MultiGridHierarchy.h"
#include "Mehrgitterverfahren/Prolongation.h"
#include "matrix/matrix.h"
#include "Mehrgitterverfahren/gmg.h"
#include "Mehrgitterverfahren/PoissonAssembler.h"

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

// int main_alt() {
//     try {
//         // --- Setup der Matrix A (Diagonaldmoninant für Konvergenz) ---
//         const size_t N_ROWS = 3;
//         const size_t N_COLS = 3;
//         const size_t MAX_ENTRIES_N = 3; 
        
//         SparseMatrix<T> A(N_ROWS, N_COLS, MAX_ENTRIES_N);

//         // Matrix A: [ 10 -1  0 ], [ -1 10 -2 ], [  0 -2 10 ]
//         A.addEntry(0, 0, 10.0);
//         A.addEntry(0, 1, -1.0);
//         A.addEntry(1, 0, -1.0);
//         A.addEntry(1, 1, 10.0);
//         A.addEntry(1, 2, -2.0);
//         A.addEntry(2, 1, -2.0);
//         A.addEntry(2, 2, 10.0);
        
//         // Vektor b: Für die bekannte Lösung x = [1, 1, 1]
//         std::vector<T> b = {9.0, 7.0, 8.0};
        
//         // Parameter
//         const int MAX_ITER = 50;
//         const T TOLERANCE = 1e-6;

//         // 1. Ausgabe der Matrix
//         printDenseMatrix(A);
        
//         std::cout << "Lösung b: [9.0, 7.0, 8.0]\n";
//         std::cout << "Gesuchte Lösung: x = [1.0, 1.0, 1.0]\n\n";

//         // --- 2. Jacobi-Verfahren ---
//         std::cout << "===========================================\n";
//         std::cout << "          STARTE JACOBI-VERFAHREN          \n";
//         std::cout << "===========================================\n";
        
//         // Wichtig: x_start wird an solveJacobi als Kopie übergeben
//         std::vector<T> x_start(N_ROWS, 0.0); 
//         std::vector<T> x_jacobi = solveJacobi<T>(A, b, x_start, MAX_ITER, TOLERANCE);
        
//         std::cout << "Jacobi FINALE LÖSUNG:\n";
//         printSolution(x_jacobi);
//         std::cout << "-------------------------------------------\n";

//         // --- 3. Gauss-Seidel-Verfahren ---
//         std::cout << "===========================================\n";
//         std::cout << "        STARTE GAUSS-SEIDEL-VERFAHREN      \n";
//         std::cout << "===========================================\n";

//         // Wichtig: Erneute Startschätzung für den fairen Vergleich
//         std::vector<T> x_start_gs(N_ROWS, 0.0);
//         std::vector<T> x_gauss_seidel = solveGaussSeidel<T>(A, b, x_start_gs, MAX_ITER, TOLERANCE);

//         std::cout << "Gauss-Seidel FINALE LÖSUNG:\n";
//         printSolution(x_gauss_seidel);
//         std::cout << "-------------------------------------------\n";


//     } catch (const std::exception& e) {
//         std::cerr << "\nFEHLER: " << e.what() << "\n";
//     }
//     return 0;
// }

// int main_Aufgabe1&2() {
//     try {
//         // --- TEST AUFGABE 1: Mehrgitterhierarchie ---
//         std::cout << "===========================================\n";
//         std::cout << "       TEST AUFGABE 1: HIERARCHIE         \n";
//         std::cout << "===========================================\n";

//         // Parameter: Einheitsquadrat (0,1)x(0,1), 3x3 Grobgitter, 1 Verfeinerung
//         // Level 0: 3x3 Knoten, Level 1: 5x5 Knoten (wegen 2n-1)
//         MultiGridHierarchy hierarchy(0.0, 1.0, 0.0, 1.0, 3, 1);

//         std::cout << "Anzahl Level: " << hierarchy.getNumLevels() << " (Erwartet: 2)\n";
//         for (int l = 0; l < hierarchy.getNumLevels(); ++l) {
//             const auto& lvl = hierarchy.getLevel(l);
//             std::cout << "Level " << l << ": " << lvl.n << "x" << lvl.n 
//                       << " Knoten, h=" << lvl.h << ", Gesamtknoten: " << lvl.numNodes << "\n";
//         }

//         // --- TEST AUFGABE 2: Prolongation ---
//         std::cout << "\n===========================================\n";
//         std::cout << "       TEST AUFGABE 2: PROLONGATION       \n";
//         std::cout << "===========================================\n";

//         Prolongation prol;
//         Matrix P(1, 1); // Wird durch assemble() resized
//         prol.assemble(P, hierarchy.getLevel(0), hierarchy.getLevel(1));

//         std::cout << "Dimension P: " << P.getRows() << "x" << P.getCols() 
//                   << " (Erwartet: 25x9)\n";

//         // Test-Vektor auf dem groben Gitter (Level 0): Alles 1.0
//         // Wenn wir eine konstante Fläche von 1.0 prolongieren (interpolieren),
//         // muss auf dem feinen Gitter (Level 1) auch überall 1.0 rauskommen.
//         std::vector<double> v_coarse(hierarchy.getLevel(0).numNodes, 1.0);
//         std::vector<double> v_fine(hierarchy.getLevel(1).numNodes, 0.0);

//         // Matrix-Vektor-Multiplikation: v_fine = P * v_coarse
//         for (size_t i = 0; i < P.getRows(); ++i) {
//             for (size_t j = 0; j < P.getCols(); ++j) {
//                 v_fine[i] += P(i, j) * v_coarse[j];
//             }
//         }

//         std::cout << "Interpolations-Check (Konstante 1.0): ";
//         bool ok = true;
//         for (double val : v_fine) {
//             if (std::abs(val - 1.0) > 1e-9) ok = false;
//         }
//         std::cout << (ok ? "ERFOLGREICH (Alle Werte sind 1.0)" : "FEHLER (Werte weichen ab)") << "\n";

//         // --- TEST HINWEIS: Restriktion ---
//         std::cout << "\n===========================================\n";
//         std::cout << "       TEST HINWEIS: RESTRIKTION          \n";
//         std::cout << "===========================================\n";
//         std::cout << "Restriktionsmatrix R = P^T (Dimension: 9x25)\n";
//         // Hier könntest du deine Transponier-Funktion testen.

//     } catch (const std::exception& e) {
//         std::cerr << "\nFEHLER: " << e.what() << "\n";
//     }
//     // ... restlicher Jacobi/GS Code
//     return 0;
// }

int main_Aufgabe3() { // MUSS int sein
    try {
        std::cout << "===========================================\n";
        std::cout << "       TEST AUFGABE 3: GMG V-ZYKLUS       \n";
        std::cout << "===========================================\n";

        // 1. Hierarchie erstellen
        MultiGridHierarchy hierarchy(0.0, 1.0, 0.0, 1.0, 3, 2);
        
        // 2. System auf dem feinsten Level vorbereiten
        PoissonAssembler assembler;
        const LevelInfo& fineLvl = hierarchy.getLevel(2);
        SparseMatrix<double> A_fine(1, 1, 5);
        std::vector<double> f_fine;
        
        assembler.assemble(A_fine, f_fine, fineLvl);

        // 3. GMG-Löser initialisieren
        GMG gmg(hierarchy);

        // Startwert u = 0
        std::vector<double> u(fineLvl.numNodes, 0.0);

        // 4. Residuum vor dem ersten V-Zyklus
        auto calculateResidual = [&](const std::vector<double>& sol) -> double {
            std::vector<double> Au = A_fine.multiply(sol);
            double resNorm = 0.0;
            for (size_t i = 0; i < f_fine.size(); ++i) {
                resNorm += (f_fine[i] - Au[i]) * (f_fine[i] - Au[i]);
            }
            return std::sqrt(resNorm);
        };

        double resStart = calculateResidual(u);
        std::cout << "Residuum am Anfang: " << resStart << "\n";

        // 5. Drei V-Zyklen durchführen
        for (int i = 1; i <= 3; ++i) {
            gmg.solve(u, f_fine);
            double resStep = calculateResidual(u);
            std::cout << "Residuum nach V-Zyklus " << i << ": " << resStep << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "GMG Test-Fehler: " << e.what() << "\n";
    }
    return 0; // Rückgabewert hinzufügen
}

// int main_Aufgabe4() {
//     std::cout << std::left << std::setw(8) << "Level" << std::setw(12) << "Knoten" 
//               << std::setw(10) << "Iter." << std::setw(15) << "Konv.rate" << std::endl;
//     std::cout << "------------------------------------------------------------" << std::endl;

//     for (int L = 1; L <= 7; ++L) {
//         MultiGridHierarchy hierarchy(0.0, 1.0, 0.0, 1.0, 3, L);
//         const LevelInfo& fineLvl = hierarchy.getLevel(L);
        
//         PoissonAssembler assembler;
//         SparseMatrix<double> A_fine;
//         std::vector<double> f_fine;
//         assembler.assemble(A_fine, f_fine, fineLvl);

//         // WICHTIG: u mit den korrekten Randwerten initialisieren
//         std::vector<double> u(fineLvl.numNodes, 0.0);
//         int n = fineLvl.n;
//         for (int i = 0; i < n; ++i) {
//             for (int j = 0; j < n; ++j) {
//                 if (i == 0 || i == n - 1 || j == 0 || j == n - 1) {
//                     double x = i * fineLvl.h;
//                     double y = j * fineLvl.h;
//                     u[i * n + j] = 1.0 - 0.5*x*x - 0.5*y*y;
//                 }
//             }
//         }

//         GMG gmg(hierarchy);
        
//         auto getNorm = [&]() {
//             std::vector<double> r_vec = A_fine.multiply(u);
//             double s = 0;
//             for(size_t i=0; i<f_fine.size(); ++i) {
//                 double diff = f_fine[i] - r_vec[i];
//                 s += diff * diff;
//             }
//             return std::sqrt(s);
//         };

//         double res = getNorm(), res_old = res, rho = 0;
//         int iter = 0;
//         // Wir lösen bis 1e-8.
//         while (res > 1e-8 && iter < 20) {
//             res_old = res;
//             gmg.solve(u, f_fine);
//             res = getNorm();
//             rho = res / res_old;
//             iter++;
//         }

//         std::cout << std::left << std::setw(8) << L << std::setw(12) << fineLvl.numNodes 
//                   << std::setw(10) << iter << std::setw(15) << std::fixed << std::setprecision(5) << rho << std::endl;
//     }
//     return 0;
// }