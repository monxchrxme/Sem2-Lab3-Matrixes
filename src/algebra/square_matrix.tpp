#pragma once

// Constructors
// We pass (size, size) to the parent constructor Matrix

template <typename T>
SquareMatrix<T>::SquareMatrix() : Matrix<T>() {}

template <typename T>
SquareMatrix<T>::SquareMatrix(int size) : Matrix<T>(size, size) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(int size, const T& fill_value)
        : Matrix<T>(size, size, fill_value) {}

template <typename T>
SquareMatrix<T>::SquareMatrix(T* items, int size)
        : Matrix<T>(items, size, size) {}


// Methods

template <typename T>
int SquareMatrix<T>::get_size() const {
    return this->get_rows();
}

template <typename T>
T SquareMatrix<T>::trace() const {
    T sum = T{};
    for (int i = 0; i < get_size(); ++i) {
        sum = sum + this->get(i, i);
    }
    return sum;
}

template <typename T>
SquareMatrix<T>* SquareMatrix<T>::transpose() const {
    int n = get_size();
    auto* result = new SquareMatrix<T>(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // swap i and j
            result->set(j, i, this->get(i, j));
        }
    }
    return result;
}

template <typename T>
bool SquareMatrix<T>::is_symmetric(double tol) const {
    int n = get_size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double diff = static_cast<double>(custom_math::abs(this->get(i, j) - this->get(j, i)));
            if (diff > tol) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
SquareMatrix<T>* SquareMatrix<T>::clone() const {
    return new SquareMatrix<T>(*this);
}

// Overriding Mathematical Operations

template <typename T>
SquareMatrix<T>* SquareMatrix<T>::add(const IMatrix<T>& other) const {
    if (this->get_rows() != other.get_rows() || this->get_cols() != other.get_cols()) {
        throw std::invalid_argument("SquareMatrix::add: Dimensions mismatch");
    }
    int n = get_size();
    auto* result = new SquareMatrix<T>(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result->set(i, j, this->get(i, j) + other.get(i, j));
        }
    }
    return result;
}

template <typename T>
SquareMatrix<T>* SquareMatrix<T>::mult(const T& value) const {
    int n = get_size();
    auto* result = new SquareMatrix<T>(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result->set(i, j, this->get(i, j) * value);
        }
    }
    return result;
}

template <typename T>
SquareMatrix<T>* SquareMatrix<T>::mult(const IMatrix<T>& other) const {
    if (this->get_cols() != other.get_rows()) {
        throw std::invalid_argument("SquareMatrix::mult: Inner dimensions must agree");
    }
    // Check: If we multiply a square matrix by a non-square one, the result will not be square
    // Therefore, this method is mathematically correct only if 'other' also has the required dimensions
    if (other.get_cols() != get_size()) {
        throw std::invalid_argument("SquareMatrix::mult: Result would not be a square matrix");
    }

    int n = get_size();
    auto* result = new SquareMatrix<T>(n, T{});
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            T temp = this->get(i, k);
            for (int j = 0; j < n; ++j) {
                result->set(i, j, result->get(i, j) + temp * other.get(k, j));
            }
        }
    }
    return result;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator+=(const IMatrix<T>& other) {
    Matrix<T>::operator+=(other); 
    return *this;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator-=(const IMatrix<T>& other) {
    Matrix<T>::operator-=(other);
    return *this;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const T& scalar) {
    Matrix<T>::operator*=(scalar);
    return *this;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator+(const IMatrix<T>& other) const {
    SquareMatrix<T> result(*this);
    result += other;
    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator-(const IMatrix<T>& other) const {
    SquareMatrix<T> result(*this);
    result -= other;
    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const T& scalar) const {
    SquareMatrix<T> result(*this);
    result *= scalar;
    return result;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const IMatrix<T>& other) const {
    if (this->get_cols() != other.get_rows()) {
        throw std::invalid_argument("SquareMatrix::operator*: Inner dimensions must agree");
    }
    if (other.get_cols() != this->get_size()) {
        throw std::invalid_argument("SquareMatrix::operator*: Result would not be a square matrix");
    }

    int n = this->get_size();
    SquareMatrix<T> result(n, T{});
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            T temp = this->get(i, k);
            for (int j = 0; j < n; ++j) {
                result.set(i, j, result.get(i, j) + temp * other.get(k, j));
            }
        }
    }
    return result;
}

template <typename T>
SquareMatrix<T> operator*(const T& scalar, const SquareMatrix<T>& m) {
    return m * scalar;
}