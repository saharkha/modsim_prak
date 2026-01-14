#ifndef PROLONGATION_H
#define PROLONGATION_H

#include "MultiGridHierarchy.h"
#include "../matrix/matrix.h" 

class Prolongation {
public:
    /**
     * Assembliert die Prolongationsmatrix P.
     * Dimension von P: (n_fine*n_fine) x (n_coarse*n_coarse)
     */
    void assemble(Matrix& P, const LevelInfo& coarse, const LevelInfo& fine);

private:
    // Hilfsfunktion zur Indexberechnung (2D -> 1D)
    int getIndex(int i, int j, int n) {
        return i * n + j;
    }
};

#endif