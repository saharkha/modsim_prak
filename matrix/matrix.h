#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>

// Wir legen den Typ auf double fest, um eine Aufteilung in .h/.cpp zu ermöglichen.
class Matrix {
private:
    std::vector<double> data;
    size_t rows;
    size_t cols;

public:
    // Konstruktor
    Matrix(size_t r, size_t c, double initial_value = 0.0);

    // --- Elementzugriff ---
    // Schreib- und Lesezugriff
    double& operator()(size_t r, size_t c);
    // Nur Lesezugriff
    const double& operator()(size_t r, size_t c) const;

    // --- Getter ---
    size_t getRows() const;
    size_t getCols() const;

    // --- Hilfsfunktion ---
    void print() const;
};

#endif // MATRIX_H