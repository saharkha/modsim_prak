#ifndef GMG_H
#define GMG_H

#include "MultiGridHierarchy.h"
#include "Prolongation.h"
#include "matrix/sparsematrix.h"
#include "matrix/matrix.h"
#include <vector>

class GMG {
public:
    GMG(MultiGridHierarchy& hierarchy);

    // Führt einen kompletten V-Zyklus auf dem feinsten Gitter aus
    void solve(std::vector<double>& u, const std::vector<double>& f);

private:
    // Der rekursive V-Zyklus gemäß Aufgabe 3
    void vCycle(int level, std::vector<double>& u, std::vector<double>& f);

    MultiGridHierarchy& m_hierarchy;
    std::vector<SparseMatrix<double>> m_A; // Systemmatrizen pro Level
    std::vector<Matrix> m_P;               // Prolongationen
    std::vector<Matrix> m_R;               // Restriktionen (R = P^T)
    
    Matrix m_baseLU; // Dichte Matrix für die LU-Zerlegung von Level 0 

    // Hilfsmethoden für LU (können auch in Matrix-Klasse oder als statische Helfer)
    void luDecompose(Matrix& A);
    std::vector<double> luSolve(const Matrix& LU, const std::vector<double>& b);
};

#endif