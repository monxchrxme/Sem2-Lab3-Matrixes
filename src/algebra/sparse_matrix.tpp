#pragma once

#include "vector.hpp"

template <typename T>
SparseMatrix<T>::SparseMatrix(int rows, int cols) : 
    m_rows(rows), 
    m_cols(cols), 
    m_data(0) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("SparseMatrix: dimensions must be positive");
    }
}

template <typename T>
int SparseMatrix<T>::get_rows() const {
    return m_rows;
}

template <typename T>
int SparseMatrix<T>::get_cols() const {
    return m_cols;
}

template <typename T>
bool SparseMatrix<T>::check_index(int i, int j) const noexcept {
    return (i >= 0 && i < m_rows && j >= 0 && j < m_cols);
}

template <typename T>
int SparseMatrix<T>::find_index(int i, int j) const {
    if (m_data.get_size() == 0) {
        return -1;
    };

    int low = 0;
    int high = m_data.get_size() - 1;

    // Binary search O(log N)
    while (low <= high) {
        int mid = low + (high - low) / 2;
        Triplet t = m_data.get(mid);

        if (t.row > i || (t.row == i && t.col > j)) {
            high = mid - 1;
        }
        else if (t.row < i || (t.row == i && t.col < j)) {
            low = mid + 1;
        }
        else {
            return mid; // found elem
        }
    }
    return -1; // not found elem
}

template <typename T>
const T& SparseMatrix<T>::get(int i, int j) const {
    if (!check_index(i, j)) {
        throw std::invalid_argument("SparseMatrix::get: index out of bounds");
    }
    
    int idx = find_index(i, j);
    if (idx != -1) {
        return m_data.get(idx).val; 
    }
    return ZERO; 
}

template <typename T>
void SparseMatrix<T>::set(int i, int j, const T& value) {
    if (!check_index(i, j)) {
        throw std::invalid_argument("SparseMatrix::set: index out of bounds");
    }

    int idx = find_index(i, j);

    if (idx != -1) {
        // Case 1: The element already exists
        if (value == T{}) {
            // Deleting: shift all elements to the right to the left
            int current_size = m_data.get_size();
            for (int k = idx; k < current_size - 1; ++k) {
                m_data.set(k, m_data.get(k + 1));
            }
            m_data.resize(current_size - 1);
        } else {
            // Simply update the value
            Triplet t = m_data.get(idx);
            t.val = value;
            m_data.set(idx, t);
        }
    } 
    else {
        // Case 2: Element not found
        if (value == T{}) {
            return; 
        }

        // Add a new element while maintaining sort order
        int current_size = m_data.get_size();
        m_data.resize(current_size + 1);
        
        int dst = current_size; // Points to the new end of the array
        
        // Shift larger elements to the right
        while (dst > 0) {
            Triplet prev = m_data.get(dst - 1);
            if (prev.row > i || (prev.row == i && prev.col > j)) {
                m_data.set(dst, prev);
                dst--;
            } else {
                break;
            }
        }
        m_data.set(dst, Triplet{i, j, value});
    }
}


template <typename T>
double SparseMatrix<T>::norm() const {
    double sum = 0.0;
    for (int k = 0; k < m_data.get_size(); ++k) {
        double abs_val = static_cast<double>(custom_math::abs(m_data.get(k).val));
        sum += abs_val * abs_val;
    }
    return std::sqrt(sum);
}

template <typename T>
SparseMatrix<T>* SparseMatrix<T>::clone() const {
    return new SparseMatrix<T>(*this);
}

template <typename T>
SparseMatrix<T>* SparseMatrix<T>::add(const IMatrix<T>& other) const {
    if (m_rows != other.get_rows() || m_cols != other.get_cols()) {
        throw std::invalid_argument("SparseMatrix::add: dimensions mismatch");
    }
    
    auto* result = new SparseMatrix<T>(m_rows, m_cols);
    
    try {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                T sum = this->get(i, j) + other.get(i, j);
                result->set(i, j, sum); 
            }
        }
    } catch (...) {
        delete result; 
        throw; 
    }
    
    return result;
}

template <typename T>
SparseMatrix<T>* SparseMatrix<T>::mult(const T& value) const {
    auto* result = new SparseMatrix<T>(m_rows, m_cols);
    if (value == T{}) return result; 

    try {
        for (int k = 0; k < m_data.get_size(); ++k) {
            Triplet t = m_data.get(k);
            result->set(t.row, t.col, t.val * value);
        }
    } catch (...) {
        delete result;
        throw;
    }
    
    return result;
}

template <typename T>
SparseMatrix<T>* SparseMatrix<T>::mult(const IMatrix<T>& other) const {
    if (m_cols != other.get_rows()) {
        throw std::invalid_argument("SparseMatrix::mult: dimensions mismatch");
    }
    
    auto* result = new SparseMatrix<T>(m_rows, other.get_cols());
    
    try {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.get_cols(); ++j) {
                T sum = T{};
                for (int k = 0; k < m_cols; ++k) {
                    sum = sum + (this->get(i, k) * other.get(k, j));
                }
                result->set(i, j, sum); 
            }
        }
    } catch (...) {
        delete result;
        throw;
    }
    
    return result;
}

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator+=(const IMatrix<T>& other) {
    if (m_rows != other.get_rows() || m_cols != other.get_cols()) {
        throw std::invalid_argument("SparseMatrix::operator+=: Dimensions mismatch");
    }

    // FAST PATH: If the other matrix is ​​also sparse, iterate only over its non-zero elements
    if (auto* sparse_other = dynamic_cast<const SparseMatrix<T>*>(&other)) {
        for (int k = 0; k < sparse_other->m_data.get_size(); ++k) {
            Triplet t = sparse_other->m_data.get(k);
            this->set(t.row, t.col, this->get(t.row, t.col) + t.val);
        }
    } 
    // SLOW PATH: Add dense matrix
    else {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                T val = other.get(i, j);
                if (val != T{}) { 
                    this->set(i, j, this->get(i, j) + val);
                }
            }
        }
    }
    return *this;
}

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator-=(const IMatrix<T>& other) {
    if (m_rows != other.get_rows() || m_cols != other.get_cols()) {
        throw std::invalid_argument("SparseMatrix::operator-=: Dimensions mismatch");
    }

    if (auto* sparse_other = dynamic_cast<const SparseMatrix<T>*>(&other)) {
        for (int k = 0; k < sparse_other->m_data.get_size(); ++k) {
            Triplet t = sparse_other->m_data.get(k);
            this->set(t.row, t.col, this->get(t.row, t.col) - t.val);
        }
    } else {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                T val = other.get(i, j);
                if (val != T{}) {
                    this->set(i, j, this->get(i, j) - val);
                }
            }
        }
    }
    return *this;
}

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator*=(const T& scalar) {
    if (scalar == T{}) {
        SparseMatrix<T> empty(m_rows, m_cols);
        *this = std::move(empty); 
        return *this;
    }

    for (int k = 0; k < m_data.get_size(); ++k) {
        Triplet t = m_data.get(k);
        t.val = t.val * scalar;
        m_data.set(k, t);
    }
    return *this;
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator+(const IMatrix<T>& other) const {
    SparseMatrix<T> result(*this);
    result += other;
    return result;
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator-(const IMatrix<T>& other) const {
    SparseMatrix<T> result(*this);
    result -= other;
    return result;
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator*(const T& scalar) const {
    SparseMatrix<T> result(*this);
    result *= scalar;
    return result;
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator*(const IMatrix<T>& other) const {
    if (m_cols != other.get_rows()) {
        throw std::invalid_argument("SparseMatrix::operator*: Inner dimensions must agree");
    }
    
    int other_cols = other.get_cols();
    SparseMatrix<T> result(m_rows, other_cols);

    for (int k = 0; k < m_data.get_size(); ++k) {
        Triplet t = m_data.get(k);
        for (int j = 0; j < other_cols; ++j) {
            T other_val = other.get(t.col, j);
            if (other_val != T{}) { 
                result.set(t.row, j, result.get(t.row, j) + t.val * other_val);
            }
        }
    }
    return result;
}

template <typename T>
Vector<T> SparseMatrix<T>::operator*(const Vector<T>& v) const {
    if (m_cols != v.get_size()) {
        throw std::invalid_argument("SparseMatrix::operator*: Dimensions mismatch");
    }
    Vector<T> result(m_rows, T{}); 

    for (int k = 0; k < m_data.get_size(); ++k) {
        Triplet t = m_data.get(k);
        result.set(t.row, result.get(t.row) + t.val * v.get(t.col));
    }
    return result;
}

template <typename T>
SparseMatrix<T> operator*(const T& scalar, const SparseMatrix<T>& m) {
    return m * scalar;
}