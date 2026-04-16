#pragma once
#include "matrix.hpp"
#include "../types/math.hpp"
#include "vector.hpp"

// auxiliary private methods (DRY)

template <typename T>
int Matrix<T>::checked_size(int rows, int cols) {
    if (rows < 0 || cols < 0) {
        throw std::invalid_argument("Matrix dimensions cannot be negative");
    }
    return rows * cols;
}

template <typename T>
bool Matrix<T>::in_bounds(int i, int j) const noexcept {
    return i >= 0 && i < m_rows && j >= 0 && j < m_cols;
}

template <typename T>
int Matrix<T>::get_index(int i, int j) const noexcept {
    return i * m_cols + j;
}

// constructors

template <typename T>
Matrix<T>::Matrix() :
    m_data(),
    m_rows(0),
    m_cols(0) {}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) :
    m_data(checked_size(rows, cols)),
    m_rows(rows),
    m_cols(cols) {}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, const T& fill_value) :
    m_data(checked_size(rows, cols)),
    m_rows(rows),
    m_cols(cols) {
    for (int i = 0; i < m_rows * m_cols; ++i) {
        m_data.set(i, fill_value);
    }
}

template <typename T>
Matrix<T>::Matrix(T* items, int rows, int cols) :
    m_data(items, checked_size(rows, cols)),
    m_rows(rows),
    m_cols(cols) {}

// Move Semantics (Moving)

// use std::move to force the DynamicArray class to
// use its move constructor rather than its copy constructor
template <typename T>
Matrix<T>::Matrix(Matrix<T>&& other) noexcept :
    m_data(std::move(other.m_data)),
    m_rows(other.m_rows),
    m_cols(other.m_cols)
{
    // Explicitly zero out the dimensions of the source matrix
    other.m_rows = 0;
    other.m_cols = 0;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& other) noexcept {
    if (this != &other) {
        m_data = std::move(other.m_data);
        m_rows = other.m_rows;
        m_cols = other.m_cols;

        // Explicitly zero out the dimensions of the source matrix
        other.m_rows = 0;
        other.m_cols = 0;
    }
    return *this;
}

// getters

template <typename T>
int Matrix<T>::get_rows() const {
    return m_rows;
}

template <typename T>
int Matrix<T>::get_cols() const {
    return m_cols;
}

template <typename T>
const T& Matrix<T>::get(int i, int j) const {
    if (!in_bounds(i, j)) {
        throw IndexOutOfRange("Matrix::get: Index out of bounds");
    }
    return m_data.get(get_index(i, j));
}

// change operations

template <typename T>
void Matrix<T>::set(int i, int j, const T& value) {
    if (!in_bounds(i, j)) {
        throw IndexOutOfRange("Matrix::set: Index out of bounds");
    }
    m_data.set(get_index(i, j), value);
}

// math operations

template <typename T>
Matrix<T>* Matrix<T>::add(const IMatrix<T>& other) const {
    if (m_rows != other.get_rows() || m_cols != other.get_cols()) {
        throw std::invalid_argument("Matrix::add: Dimensions mismatch");
    }
    auto* result = new Matrix<T>(m_rows, m_cols);
    for (int i = 0; i < m_rows; ++i) {
        for (int j = 0; j < m_cols; ++j) {
            result->set(i, j, this->get(i, j) + other.get(i, j));
        }
    }
    return result;
}

template <typename T>
Matrix<T>* Matrix<T>::mult_scalar(const T& value) const {
    auto* result = new Matrix<T>(m_rows, m_cols);
    for (int i = 0; i < m_rows * m_cols; ++i) {
        result->m_data.set(i, m_data.get(i) * value);
    }
    return result;
}

template <typename T>
Matrix<T>* Matrix<T>::mult(const IMatrix<T>& other) const {
    if (m_cols != other.get_rows()) {
        throw std::invalid_argument("Matrix::mult: Inner dimensions must agree");
    }
    int other_cols = other.get_cols();
    auto* result = new Matrix<T>(m_rows, other_cols, T{}); // initializing with zeroes (T{})

    for (int i = 0; i < m_rows; ++i) {
        for (int k = 0; k < m_cols; ++k) {
            T temp = this->get(i, k);
            for (int j = 0; j < other_cols; ++j) {
                result->set(i, j, result->get(i, j) + temp * other.get(k, j));
            }
        }
    }
    return result;
}

template <typename T>
Vector<T> Matrix<T>::operator*(const Vector<T>& v) const {
    if (m_cols != v.get_size()) {
        throw std::invalid_argument("Matrix-Vector multiplication: Dimensions mismatch");
    }
    Vector<T> result(m_rows, T{}); // Create a vector of size m_rows, filled with zeros

    for (int i = 0; i < m_rows; ++i) {
        T sum = T{};
        for (int j = 0; j < m_cols; ++j) {
            sum = sum + this->get(i, j) * v.get(j);
        }
        result.set(i, sum);
    }
    return result;
}

template <typename T>
double Matrix<T>::norm() const {
    double sum = 0;
    for (int i = 0; i < m_rows * m_cols; ++i) {
        double a = static_cast<double>(custom_math::abs(m_data.get(i)));
        sum += a * a;
    }
    return std::sqrt(sum);
}

template <typename T>
Matrix<T>* Matrix<T>::clone() const {
    return new Matrix<T>(*this);
}

// specific methods for linear algebra

template <typename T>
void Matrix<T>::swap_rows(int row1, int row2) {
    if (row1 == row2) return;
    if (row1 < 0 || row1 >= m_rows || row2 < 0 || row2 >= m_rows) {
        throw IndexOutOfRange("Matrix::swap_rows: Index out of bounds");
    }
    for (int j = 0; j < m_cols; ++j) {
        T temp = m_data.get(get_index(row1, j));
        m_data.set(get_index(row1, j), m_data.get(get_index(row2, j)));
        m_data.set(get_index(row2, j), temp);
    }
}

template <typename T>
void Matrix<T>::row_operation(int target_row, int source_row, const T& factor) {
    if (target_row < 0 || target_row >= m_rows || source_row < 0 || source_row >= m_rows) {
        throw IndexOutOfRange("Matrix::row_operation: Index out of bounds");
    }
    for (int j = 0; j < m_cols; ++j) {
        T current = m_data.get(get_index(target_row, j));
        T source_val = m_data.get(get_index(source_row, j));
        m_data.set(get_index(target_row, j), current - factor * source_val);
    }
}

template <typename T>
void Matrix<T>::scale_row(int row, const T& factor) {
    if (row < 0 || row >= m_rows) {
        throw IndexOutOfRange("scale_row: Index out of bounds");
    }
    for (int j = 0; j < m_cols; ++j) {
        m_data.set(get_index(row, j), m_data.get(get_index(row, j)) * factor);
    }
}