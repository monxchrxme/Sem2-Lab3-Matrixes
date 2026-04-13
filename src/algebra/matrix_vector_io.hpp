#pragma once

#include <iostream>
#include <iomanip>
#include "imatrix.hpp"
#include "vector.hpp"

template <typename T>
std::ostream& operator<<(std::ostream& os, const IMatrix<T>& m) {
    for (int row = 0; row < m.get_rows(); ++row) {
        os << "( ";
        for (int col = 0; col < m.get_cols(); ++col) {
            os << std::setw(8) << m.get(row, col);

            if (col < m.get_cols() - 1) {
                os << ", ";
            }
        }
        os << " )\n";
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
    os << "[ ";
    for (int i = 0; i < v.get_size(); ++i) {
        os << v.get(i);
        if (i < v.get_size() - 1) {
            os << ", ";
        }
    }
    os << " ]";
    return os;
}