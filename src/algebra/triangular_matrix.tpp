#pragma once
#include "triangular_matrix.hpp"
#include "../types/math.hpp"

// Constructors

template <typename T>
TriangularMatrix<T>::TriangularMatrix(TriangleType type) :
    m_data(),
    m_size(0),
    m_type(type) {}

template <typename T>
TriangularMatrix<T>::TriangularMatrix(int size, TriangleType type) :
    m_data(size * (size + 1) / 2),
    m_size(size),
    m_type(type) {
    if (size < 0) {
        throw std::invalid_argument("Size cannot be negative");
    }
}

template <typename T>
TriangularMatrix<T>::TriangularMatrix(int size, const T& fill_value, TriangleType type) :
    m_data(size * (size + 1) / 2),
    m_size(size),
    m_type(type) {
    if (size < 0) {
        throw std::invalid_argument("Size cannot be negative");
    }
    for (int i = 0; i < m_data.get_size(); ++i) {
        m_data.set(i, fill_value);
    }
}

// Movement (As on Diagonal)

template <typename T>
TriangularMatrix<T>::TriangularMatrix(TriangularMatrix<T>&& other) noexcept :
    m_data(std::move(other.m_data)),
    m_size(other.m_size),
    m_type(other.m_type) {
    other.m_size = 0;
}

template <typename T>
TriangularMatrix<T>& TriangularMatrix<T>::operator=(TriangularMatrix<T>&& other) noexcept {
    if (this != &other) {
        m_data = std::move(other.m_data);
        m_size = other.m_size;
        m_type = other.m_type;
        other.m_size = 0;
    }
    return *this;
}

// getters

template <typename T>
int TriangularMatrix<T>::get_rows() const {
    return m_size;
}

template <typename T>
int TriangularMatrix<T>::get_cols() const {
    return m_size;
}

template <typename T>
int TriangularMatrix<T>::get_size() const {
    return m_size;
}

template <typename T>
TriangleType TriangularMatrix<T>::get_type() const {
    return m_type;
}

// Access to elements

template <typename T>
const T& TriangularMatrix<T>::get(int i, int j) const {
    if (!in_bounds(i, j)) {
        throw IndexOutOfRange("TriangularMatrix::get: Index out of bounds");
    }

    if (!in_structure(i, j)) {
        return ZERO; // If in the zero zone, return our fast constant
    }
    return m_data.get(flat_index(i, j));
}

template <typename T>
void TriangularMatrix<T>::set(int i, int j, const T& value) {
    if (!in_bounds(i, j)) {
        throw IndexOutOfRange("TriangularMatrix::set: Index out of bounds");
    }

    if (!in_structure(i, j)) {
        if (value != T{}) {
            throw std::logic_error("Cannot set non-zero value outside the triangle structure");
        }
        return; // Ignore writing zero to the zero zone
    }
    m_data.set(flat_index(i, j), value);
}

// Math

template <typename T>
SquareMatrix<T>* TriangularMatrix<T>::add(const IMatrix<T>& other) const {
    if (m_size != other.get_rows() || m_size != other.get_cols()) {
        throw std::invalid_argument("TriangularMatrix::add: Dimensions mismatch");
    }
    auto* result = new SquareMatrix<T>(m_size);
    for (int i = 0; i < m_size; ++i) {
        for (int j = 0; j < m_size; ++j) {
            result->set(i, j, this->get(i, j) + other.get(i, j));
        }
    }
    return result;
}

template <typename T>
TriangularMatrix<T>* TriangularMatrix<T>::mult_scalar(const T& value) const {
    auto* result = new TriangularMatrix<T>(m_size, m_type);
    for (int i = 0; i < m_data.get_size(); ++i) {
        result->m_data.set(i, m_data.get(i) * value);
    }
    return result;
}

template <typename T>
SquareMatrix<T>* TriangularMatrix<T>::mult(const IMatrix<T>& other) const {
    if (m_size != other.get_rows() || m_size != other.get_cols()) {
        throw std::invalid_argument("TriangularMatrix::mult: Dimensions mismatch");
    }

    auto* result = new SquareMatrix<T>(m_size, T{});

    for (int i = 0; i < m_size; ++i) {
        for (int k = 0; k < m_size; ++k) {
            T temp = this->get(i, k);
            if (temp == T{}) continue; // Optimization: skip zero elements

            for (int j = 0; j < m_size; ++j) {
                result->set(i, j, result->get(i, j) + temp * other.get(k, j));
            }
        }
    }
    return result;
}

template <typename T>
double TriangularMatrix<T>::norm() const {
    double sum = 0;
    // Optimization: Iterate only over significant elements
    for (int i = 0; i < m_data.get_size(); ++i) {
        double a = static_cast<double>(custom_math::abs(m_data.get(i)));
        sum += a * a;
    }
    return std::sqrt(sum);
}

template <typename T>
TriangularMatrix<T>* TriangularMatrix<T>::clone() const {
    return new TriangularMatrix<T>(*this);
}

// Specific Metrics (O(N))

template <typename T>
T TriangularMatrix<T>::trace() const {
    T sum = T{};
    for (int i = 0; i < m_size; ++i) {
        sum = sum + this->get(i, i);
    }
    return sum;
}

template <typename T>
T TriangularMatrix<T>::determinant() const {
    if (m_size == 0) return T{};
    T det = T(1);
    for (int i = 0; i < m_size; ++i) {
        det = det * this->get(i, i);
    }
    return det;
}

template <typename T>
bool TriangularMatrix<T>::in_bounds(int i, int j) const noexcept {
    return i >= 0 && i < m_size && j >= 0 && j < m_size;
}

template <typename T>
bool TriangularMatrix<T>::in_structure(int i, int j) const noexcept {
    if (m_type == TriangleType::Lower) {
        return i >= j; // Significant elements below and on the diagonal
    } else {
        return i <= j; // Significant elements above and on the diagonal
    }
}

template <typename T>
int TriangularMatrix<T>::flat_index(int i, int j) const noexcept {
    // implementing our unified formula: x*(x+1)/2 + y
    if (m_type == TriangleType::Lower) {
        return i * (i + 1) / 2 + j;
    } else {
        return j * (j + 1) / 2 + i; // just swap i and j
    }
}