#pragma once
#include "diagonal_matrix.hpp"

// Constructors

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix() : m_data(), m_size(0) {}

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(int size) : m_data(size), m_size(size) {
    if (size < 0) throw std::invalid_argument("DiagonalMatrix: size cannot be negative");
}

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(int size, const T& fill_value) :
    m_data(size),
    m_size(size) {
    if (size < 0) throw std::invalid_argument("DiagonalMatrix: size cannot be negative");
    for (int i = 0; i < size; ++i) {
        m_data.set(i, fill_value);
    }
}

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(T* diag_items, int size) :
        m_data(diag_items, size),
        m_size(size) {
    if (size < 0) throw std::invalid_argument("DiagonalMatrix: size cannot be negative");
}

// move-semantics

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(DiagonalMatrix<T>&& other) noexcept :
        m_data(std::move(other.m_data)),
        m_size(other.m_size) {
    other.m_size = 0;
}

template <typename T>
DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(DiagonalMatrix<T>&& other) noexcept {
    if (this != &other) {
        m_data = std::move(other.m_data);
        m_size = other.m_size;

        other.m_size = 0;
    }
    return *this;
}

// getters

template <typename T>
int DiagonalMatrix<T>::get_rows() const { return m_size; }

template <typename T>
int DiagonalMatrix<T>::get_cols() const { return m_size; }

template <typename T>
int DiagonalMatrix<T>::get_size() const { return m_size; }

template <typename T>
const T& DiagonalMatrix<T>::get(int i, int j) const {
    if (i < 0 || i >= m_size || j < 0 || j >= m_size) {
        throw IndexOutOfRange("DiagonalMatrix::get: Index out of bounds");
    }
    if (i == j) {
        return m_data.get(i);
    }
    return ZERO;
}

template <typename T>
void DiagonalMatrix<T>::set(int i, int j, const T& value) {
    if (i < 0 || i >= m_size || j < 0 || j >= m_size) {
        throw IndexOutOfRange("DiagonalMatrix::set: Index out of bounds");
    }
    if (i == j) {
        m_data.set(i, value);
    } else {
        if (value != T{}) {
            throw std::logic_error("Cannot set non-zero value outside the main diagonal of DiagonalMatrix");
        }
    }
}

// Optimized Mathematics

template <typename T>
SquareMatrix<T>* DiagonalMatrix<T>::add(const IMatrix<T>& other) const {
    if (m_size != other.get_rows() || m_size != other.get_cols()) {
        throw std::invalid_argument("DiagonalMatrix::add: Dimensions mismatch");
    }
    auto* result = new SquareMatrix<T>(m_size, m_size);
    for (int i = 0; i < m_size; ++i) {
        for (int j = 0; j < m_size; ++j) {
            result->set(i, j, this->get(i, j) + other.get(i, j));
        }
    }
    return result;
}

template <typename T>
DiagonalMatrix<T>* DiagonalMatrix<T>::mult_scalar(const T& value) const {
    auto* result = new DiagonalMatrix<T>(m_size);
    for (int i = 0; i < m_size; ++i) {
        result->m_data.set(i, m_data.get(i) * value);
    }
    return result;
}

template <typename T>
SquareMatrix<T>* DiagonalMatrix<T>::mult(const IMatrix<T>& other) const {
    if (m_size != other.get_rows()) {
        throw std::invalid_argument("DiagonalMatrix::mult: Inner dimensions must agree");
    }
    int other_cols = other.get_cols();
    auto* result = new SquareMatrix<T>(m_size, other_cols);

    for (int i = 0; i < m_size; ++i) {
        T diag_val = m_data.get(i);
        for (int j = 0; j < other_cols; ++j) {
            result->set(i, j, diag_val * other.get(i, j));
        }
    }
    return result;
}

template <typename T>
double DiagonalMatrix<T>::norm() const {
    double sum = 0;
    for (int i = 0; i < m_size; ++i) {
        double a = static_cast<double>(custom_math::abs(m_data.get(i)));
        sum += a * a;
    }
    return std::sqrt(sum);
}

template <typename T>
DiagonalMatrix<T>* DiagonalMatrix<T>::clone() const {
    return new DiagonalMatrix<T>(*this);
}

template <typename T>
T DiagonalMatrix<T>::trace() const {
    T sum = T{};
    for (int i = 0; i < m_size; ++i) {
        sum = sum + m_data.get(i);
    }
    return sum;
}

template <typename T>
T DiagonalMatrix<T>::determinant() const {
    if (m_size == 0) return T{};
    T det = T(1);
    for (int i = 0; i < m_size; ++i) {
        det = det * m_data.get(i);
    }
    return det;
}