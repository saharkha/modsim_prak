#include "MultiGridHierarchy.h"

MultiGridHierarchy::MultiGridHierarchy(double xMin, double xMax, double yMin, double yMax, 
                                       int coarseN, int numRefinements)
    : m_xMin(xMin), m_xMax(xMax), m_yMin(yMin), m_yMax(yMax) 
{
    // 1. Grobes Gitter (Level 0) erstellen 
    LevelInfo coarse;
    coarse.n = coarseN;
    coarse.h = (xMax - xMin) / (coarseN - 1);
    coarse.numNodes = coarseN * coarseN;
    m_levels.push_back(coarse);

    // 2. Gewünschte Anzahl an Verfeinerungen durchführen 
    for (int i = 0; i < numRefinements; ++i) {
        refine();
    }
}

void MultiGridHierarchy::refine() {
    const LevelInfo& fine = m_levels.back();
    
    LevelInfo next;
    // Formel: Wenn wir h halbieren, verdoppeln wir die Intervalle.
    // Knoten = Intervalle + 1 -> $n_{neu} = 2(n-1) + 1 = 2n - 1$
    next.n = 2 * fine.n - 1;
    next.h = fine.h / 2.0;
    next.numNodes = next.n * next.n;
    
    m_levels.push_back(next);
}

int MultiGridHierarchy::getNumLevels() const {
    return static_cast<int>(m_levels.size());
}

const LevelInfo& MultiGridHierarchy::getLevel(int lvl) const {
    return m_levels.at(lvl);
}