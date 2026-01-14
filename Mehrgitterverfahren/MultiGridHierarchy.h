#ifndef MULTIGRID_HIERARCHY_H
#define MULTIGRID_HIERARCHY_H

#include <vector>

// Speichert Informationen für eine einzelne Gitterebene
struct LevelInfo {
    int n;          // Anzahl der Knoten pro Dimension (n x n Gitter)
    double h;       // Gitterweite (Abstand zwischen zwei Knoten)
    int numNodes;   // Gesamtzahl der Knoten im 2D-Gitter (n * n)
};

class MultiGridHierarchy {
public:
    /**
     * @param xMin, xMax: Räumliche Ausdehnung in x-Richtung
     * @param yMin, yMax: Räumliche Ausdehnung in y-Richtung
     * @param coarseN: Anzahl der Knoten pro Dimension auf dem gröbsten Gitter (Level 0)
     * @param numRefinements: Wie oft das Gitter verfeinert werden soll
     */
    MultiGridHierarchy(double xMin, double xMax, double yMin, double yMax, 
                       int coarseN, int numRefinements);

    // Hilfsmethoden
    int getNumLevels() const;
    const LevelInfo& getLevel(int lvl) const;

private:
    // Erzeugt ein neues, feineres Gitter basierend auf dem aktuell feinsten
    void refine();

    std::vector<LevelInfo> m_levels;
    double m_xMin, m_xMax, m_yMin, m_yMax;
};

#endif