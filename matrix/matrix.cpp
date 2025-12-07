#include "matrix.h"

// Konstruktor
Matrix::Matrix(size_t r, size_t c, double initial_value)
    : rows(r), cols(c) {

    if (rows == 0 || cols == 0) {
        throw std::invalid_argument("Matrix dimensions must be positive.");
    }
    // Reserviert Speicher für rows * cols Elemente und initialisiert sie
    data.resize(rows * cols, initial_value);
}

// --- Elementzugriff (Schreib- und Lesezugriff) ---
double& Matrix::operator()(size_t r, size_t c) {
    // Optional: Bounds-Checks hinzufügen
    if (r >= rows || c >= cols) {
        throw std::out_of_range("Matrix index out of bounds.");
    }
    // Row-Major Order: Element bei (r, c) liegt bei r * cols + c
    return data[r * cols + c];
}

// --- Elementzugriff (Nur Lesezugriff) ---
const double& Matrix::operator()(size_t r, size_t c) const {
    if (r >= rows || c >= cols) {
        throw std::out_of_range("Matrix index out of bounds.");
    }
    return data[r * cols + c];
}

// --- Getter ---
size_t Matrix::getRows() const {
    return rows;
}

size_t Matrix::getCols() const {
    return cols;
}

// --- Hilfsfunktion zur Ausgabe der Matrix ---
void Matrix::print() const {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            // Wir nutzen den konstanten operator()
            std::cout << (*this)(i, j) << "\t";
        }
        std::cout << "\n";
    }
}