#include "Prolongation.h"

void Prolongation::assemble(Matrix& P, const LevelInfo& coarse, const LevelInfo& fine) {
    // Matrixgröße anpassen: Zeilen = feine Knoten, Spalten = grobe Knoten
    P.resize(fine.numNodes, coarse.numNodes, 0.0);

    int nC = coarse.n; // Knoten pro Dimension (grob)
    int nF = fine.n;   // Knoten pro Dimension (fein)

    for (int i = 0; i < nC; ++i) {
        for (int j = 0; j < nC; ++j) {
            // Grober Index
            int cIdx = getIndex(i, j, nC);
            
            // Entsprechender Index auf dem feinen Gitter
            int fi = 2 * i;
            int fj = 2 * j;

            // 1. Der Punkt selbst (Gewicht 1.0)
            P(getIndex(fi, fj, nF), cIdx) = 1.0;

            // 2. Horizontale Nachbarn auf dem feinen Gitter (Gewicht 0.5)
            if (fi + 1 < nF) {
                P(getIndex(fi + 1, fj, nF), cIdx) = 0.5;
            }
            if (fi - 1 >= 0) {
                P(getIndex(fi - 1, fj, nF), cIdx) = 0.5;
            }

            // 3. Vertikale Nachbarn auf dem feinen Gitter (Gewicht 0.5)
            if (fj + 1 < nF) {
                P(getIndex(fi, fj + 1, nF), cIdx) = 0.5;
            }
            if (fj - 1 >= 0) {
                P(getIndex(fi, fj - 1, nF), cIdx) = 0.5;
            }

            // 4. Diagonale Nachbarn (Zentrum einer groben Zelle) (Gewicht 0.25)
            if (fi + 1 < nF && fj + 1 < nF) P(getIndex(fi + 1, fj + 1, nF), cIdx) = 0.25;
            if (fi + 1 < nF && fj - 1 >= 0) P(getIndex(fi + 1, fj - 1, nF), cIdx) = 0.25;
            if (fi - 1 >= 0 && fj + 1 < nF) P(getIndex(fi - 1, fj + 1, nF), cIdx) = 0.25;
            if (fi - 1 >= 0 && fj - 1 >= 0) P(getIndex(fi - 1, fj - 1, nF), cIdx) = 0.25;
        }
    }
}