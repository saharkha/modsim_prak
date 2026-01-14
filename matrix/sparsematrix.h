#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iterator> // Für std::forward_iterator_tag

template <typename T>
class SparseMatrix {
public:

    SparseMatrix() : rows(0), cols(0), max_entries_per_row(0) {}

    // Konstante für den Füllwert in colInds
    static constexpr size_t FILL_INDEX = (size_t)-1;
    
    // --- 1. Definition der Elemente und Iteratoren ---
    
    // Struktur, die ein Nicht-Null-Element repräsentiert (Wert und Spaltenindex)
    struct Entry {
        T value;
        size_t col_index;
    };

    class EntryIterator {
    private:
        const T* current_value_ptr;
        const size_t* current_col_ptr;

    public:
        // Standard-Iterator-Typen für STL-Kompatibilität
        using iterator_category = std::forward_iterator_tag;
        using value_type = Entry;
        using difference_type = std::ptrdiff_t;
        using pointer = const Entry*;
        using reference = const Entry&;

        // Konstruktor
        EntryIterator(const T* val_ptr, const size_t* col_ptr)
            : current_value_ptr(val_ptr), current_col_ptr(col_ptr) {}

        // Dereferenzierungsoperator: Liefert den aktuellen Eintrag
        value_type operator*() const {
            return {*current_value_ptr, *current_col_ptr};
        }

        // Prä-Inkrement-Operator: Geht zum nächsten Eintrag
        EntryIterator& operator++() {
            ++current_value_ptr;
            ++current_col_ptr;
            return *this;
        }

        // Post-Inkrement-Operator
        EntryIterator operator++(int) {
            EntryIterator temp = *this;
            ++(*this);
            return temp;
        }

        // Vergleichsoperatoren
        bool operator==(const EntryIterator& other) const {
            return current_value_ptr == other.current_value_ptr;
        }

        bool operator!=(const EntryIterator& other) const {
            return !(*this == other);
        }
    };
    
private:
    // Dimensionen der Matrix
    size_t rows;
    size_t cols;
    
    // Maximale Anzahl speicherbarer Nicht-Null-Einträge pro Zeile (N oder A)
    size_t max_entries_per_row; 

    // Speichervektoren (Zeilenweise flachgelegt)
    std::vector<size_t> colInds; 
    std::vector<T> values;

public:
    // --- 2. Konstruktor und Getter ---

    SparseMatrix(size_t r, size_t c, size_t n_max) 
        : rows(r), cols(c), max_entries_per_row(n_max) {
        
        size_t total_size = rows * max_entries_per_row;
        
        colInds.resize(total_size, FILL_INDEX); 
        values.resize(total_size, 0.0);
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    size_t getMaxEntriesPerRow() const { return max_entries_per_row; }

    // --- 3. Kernfunktionalität (Implementierung der Methoden) ---

    // Element-Zugriff A[r][c] (ACHTUNG: Langsam, erfordert Suche)
    T operator()(size_t r, size_t c) const {
        if (r >= rows || c >= cols) {
            throw std::out_of_range("Index außerhalb der Matrixgrenzen.");
        }
        
        size_t start_index = r * max_entries_per_row;
        
        for (size_t i = 0; i < max_entries_per_row; ++i) {
            size_t current_index = start_index + i;
            size_t current_col = colInds[current_index];

            if (current_col == c) {
                return values[current_index];
            }
            
            // Suche abbrechen, wenn sortiert nicht gefunden oder Füll-Index erreicht
            if (current_col == FILL_INDEX || current_col > c) {
                 return 0.0;
            }
        }
        return 0.0;
    }

    // --- Neue Resize-Methode für die Gitterhierarchie ---
    /**
     * @param r      Neue Zeilenanzahl
     * @param c      Neue Spaltenanzahl
     * @param n_max  Maximale Einträge pro Zeile (für Poisson in 2D ist das 5)
     */
    void resize(size_t r, size_t c, size_t n_max) {
        rows = r;
        cols = c;
        max_entries_per_row = n_max;
        
        size_t total_size = rows * max_entries_per_row;
        
        // Bestehende Daten löschen und mit Initialwerten füllen
        colInds.assign(total_size, FILL_INDEX); 
        values.assign(total_size, static_cast<T>(0.0));
    }
    
    // Fügt einen Eintrag hinzu und hält die Sortierung bei
    void addEntry(size_t r, size_t c, T val) {
        if (val == (T)0.0) return;
        if (r >= rows || c >= cols) return; 

        size_t start_index = r * max_entries_per_row;
        size_t insert_pos = max_entries_per_row; 
        
        // 1. Finde den Platz zum Einfügen/Aktualisieren
        for (size_t i = 0; i < max_entries_per_row; ++i) {
            size_t current_flat_index = start_index + i;
            size_t current_col = colInds[current_flat_index];

            if (current_col == c) {
                values[current_flat_index] = val;
                return;
            }
            
            if (current_col == FILL_INDEX || current_col > c) {
                insert_pos = i; 
                break;
            }
        }

        // 2. Element einfügen (Verschiebung notwendig)
        if (insert_pos < max_entries_per_row) {
            size_t insert_flat_index = start_index + insert_pos;

            // Verschiebe alle nachfolgenden Elemente (außer das letzte, das verloren geht)
            for (size_t i = max_entries_per_row - 1; i > insert_pos; --i) {
                size_t dest_idx = start_index + i;
                size_t src_idx = start_index + i - 1;
                
                colInds[dest_idx] = colInds[src_idx];
                values[dest_idx] = values[src_idx];
            }
            
            // Füge das neue Element ein
            colInds[insert_flat_index] = c;
            values[insert_flat_index] = val;
        }
    }

    // Matrix-Vektor-Multiplikation (A * x)
    std::vector<T> multiply(const std::vector<T>& x) const {
        if (cols != x.size()) {
            throw std::runtime_error("Matrix column count does not match vector size.");
        }

        std::vector<T> y(rows, 0.0);

        for (size_t r = 0; r < rows; ++r) {
            size_t start_index = r * max_entries_per_row;
            T row_sum = 0.0;

            for (size_t i = 0; i < max_entries_per_row; ++i) {
                size_t current_flat_index = start_index + i;
                size_t col_index = colInds[current_flat_index];
                
                if (col_index == FILL_INDEX) {
                    break; 
                }
                
                // Bounds Check for x is implicitly handled by the initial check, 
                // but good practice if col_index was potentially larger than cols
                if (col_index >= x.size()) {
                    // Sollte nicht passieren, wenn addEntry korrekt verwendet wurde
                    continue; 
                }

                T value = values[current_flat_index];
                row_sum += value * x[col_index];
            }
            y[r] = row_sum;
        }
        return y;
    }
    
    // --- 4. Iterator-Zugriffsmethoden ---

    // Gibt den Start-Iterator für Zeile 'r' zurück
    EntryIterator row_begin(size_t r) const {
        if (r >= rows) throw std::out_of_range("Row index out of bounds.");
        size_t start_index = r * max_entries_per_row;
        
        const T* val_ptr = values.data() + start_index;
        const size_t* col_ptr = colInds.data() + start_index;

        return EntryIterator(val_ptr, col_ptr);
    }

    // Gibt den End-Iterator für Zeile 'r' zurück (zeigt hinter das letzte gespeicherte Element)
    EntryIterator row_end(size_t r) const {
        if (r >= rows) throw std::out_of_range("Row index out of bounds.");
        size_t end_index = (r + 1) * max_entries_per_row;

        const T* val_ptr = values.data() + end_index;
        const size_t* col_ptr = colInds.data() + end_index;
        
        return EntryIterator(val_ptr, col_ptr);
    }
};

// Am Ende der sparsematrix.h, außerhalb der class SparseMatrix { ... };
template <typename T>
constexpr size_t SparseMatrix<T>::FILL_INDEX;

#endif // SPARSEMATRIX_H