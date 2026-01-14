#include "MultiGridHierarchy.h"
#include "matrix/sparsematrix.h"
#include <vector>
#include <cmath>

class PoissonAssembler {
public:
    // Assembliert A und f für ein bestimmtes Level
    void assemble(SparseMatrix<double>& A, std::vector<double>& f, const LevelInfo& lvl) {
        int n = lvl.n;
        double h = lvl.h;
        double h2_inv = 1.0 / (h * h);

        A.resize(lvl.numNodes, lvl.numNodes, 5); // Max 5 Einträge pro Zeile
        f.assign(lvl.numNodes, 2.0); // f(x,y) = 2 gemäß Aufgabe 4

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int row = i * n + j;
                double x = i * h; // Annahme: Einheitsquadrat startet bei 0.0
                double y = j * h;

                // Randknoten prüfen (Dirichlet-Randbedingung)
                if (i == 0 || i == n - 1 || j == 0 || j == n - 1) {
                    A.addEntry(row, row, 1.0);
                    f[row] = 1.0 - 0.5 * x * x - 0.5 * y * y; // u = g auf Rand 
                } 
                else {
                    // Innere Punkte: 5-Punkt-Stern
                    A.addEntry(row, row, 4.0 * h2_inv);
                    A.addEntry(row, (i - 1) * n + j, -1.0 * h2_inv);
                    A.addEntry(row, (i + 1) * n + j, -1.0 * h2_inv);
                    A.addEntry(row, i * n + (j - 1), -1.0 * h2_inv);
                    A.addEntry(row, i * n + (j + 1), -1.0 * h2_inv);
                }
            }
        }
    }
};