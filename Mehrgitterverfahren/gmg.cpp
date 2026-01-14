#include "gmg.h"
#include "PoissonAssembler.h"
#include "verfahren/GaussSeidel.h"
#include <iostream>
#include <cmath>

GMG::GMG(MultiGridHierarchy& hierarchy) : m_hierarchy(hierarchy), m_baseLU(1, 1) {
    int numLevels = m_hierarchy.getNumLevels();
    m_A.resize(numLevels);
    m_P.resize(numLevels - 1, Matrix(1, 1));
    m_R.resize(numLevels - 1, Matrix(1, 1));
    
    PoissonAssembler assembler;
    Prolongation prol;

    for (int l = 0; l < numLevels; ++l) {
        std::vector<double> dummy_f;
        assembler.assemble(m_A[l], dummy_f, m_hierarchy.getLevel(l));
        if (l < numLevels - 1) {
            prol.assemble(m_P[l], m_hierarchy.getLevel(l), m_hierarchy.getLevel(l+1));
            m_R[l].resize(m_P[l].getCols(), m_P[l].getRows());
            for (size_t r = 0; r < m_P[l].getRows(); ++r)
                for (size_t c = 0; c < m_P[l].getCols(); ++c)
                    m_R[l](c, r) = m_P[l](r, c);
        }
    }

    size_t n0 = m_A[0].getRows();
    m_baseLU.resize(n0, n0);
    for (size_t i = 0; i < n0; ++i)
        for (size_t j = 0; j < n0; ++j)
            m_baseLU(i, j) = m_A[0](i, j);
    luDecompose(m_baseLU);
}

void GMG::vCycle(int level, std::vector<double>& u, std::vector<double>& f) {
    if (level == 0) {
        u = luSolve(m_baseLU, f);
        return;
    }

    // 1. Vorglättung
    u = solveGaussSeidel(m_A[level], f, u, 2, 1e-15);

    // 2. Residuum
    std::vector<double> Au = m_A[level].multiply(u);
    std::vector<double> d(f.size());
    for (size_t i = 0; i < f.size(); ++i) d[i] = f[i] - Au[i];

    // 3. Restriktion
    size_t coarseNodes = m_hierarchy.getLevel(level - 1).numNodes;
    std::vector<double> f_coarse(coarseNodes, 0.0);
    const Matrix& R = m_R[level-1];
    for (size_t i = 0; i < R.getRows(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < R.getCols(); ++j) sum += R(i, j) * d[j];
        f_coarse[i] = sum * 0.0625; // Skalierung 1/16
    }

    // 4. Grobgitterkorrektur (Rekursion)
    std::vector<double> c_coarse(coarseNodes, 0.0);
    vCycle(level - 1, c_coarse, f_coarse);

    // 5. Prolongation
    std::vector<double> c_fine(u.size(), 0.0);
    const Matrix& P = m_P[level-1];
    for (size_t i = 0; i < P.getRows(); ++i)
        for (size_t j = 0; j < P.getCols(); ++j) c_fine[i] += P(i, j) * c_coarse[j];

    // 6. Update
    for (size_t i = 0; i < u.size(); ++i) u[i] += c_fine[i];

    // 7. Nachglättung
    u = solveGaussSeidel(m_A[level], f, u, 2, 1e-15);
}

void GMG::solve(std::vector<double>& u, const std::vector<double>& f) {
    vCycle(m_hierarchy.getNumLevels() - 1, u, const_cast<std::vector<double>&>(f));
}

// LU-Methoden
void GMG::luDecompose(Matrix& A) {
    size_t n = A.getRows();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            A(j, i) /= A(i, i);
            for (size_t k = i + 1; k < n; k++) A(j, k) -= A(j, i) * A(i, k);
        }
    }
}

std::vector<double> GMG::luSolve(const Matrix& LU, const std::vector<double>& b) {
    size_t n = LU.getRows();
    std::vector<double> x(n), y(n);
    for (size_t i = 0; i < n; i++) {
        y[i] = b[i];
        for (size_t j = 0; j < i; j++) y[i] -= LU(i, j) * y[j];
    }
    for (int i = (int)n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (size_t j = (size_t)i + 1; j < n; j++) x[i] -= LU(i, j) * x[j];
        x[i] /= LU(i, i);
    }
    return x;
}