#include "sparse_matrix.hpp"


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
bool SparseMatrix<T>::in_bounds(int i, int j) const noexcept {
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
    if (!in_bounds(i, j)) {
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
    if (!in_bounds(i, j)) {
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
SparseMatrix<T>* SparseMatrix<T>::mult_scalar(const T& value) const {
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