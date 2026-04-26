#pragma once

#include "matrix.hpp"
#include "../types/math.hpp"

template <typename T>
class SquareMatrix : public Matrix<T> {
private: 
    [[nodiscard]] int get_rows() const override { return Matrix<T>::get_rows(); }
    [[nodiscard]] int get_cols() const override { return Matrix<T>::get_cols(); }
    
public:
    // Constructors (strictly 1 size, to guarantee a square)
    SquareMatrix();
    explicit SquareMatrix(int size);
    SquareMatrix(int size, const T& fill_value);
    SquareMatrix(T* items, int size);

    // Rule of Zero
    ~SquareMatrix() override = default;
    SquareMatrix(const SquareMatrix<T>& other) = default;
    SquareMatrix(SquareMatrix<T>&& other) noexcept = default;
    SquareMatrix<T>& operator=(const SquareMatrix<T>& other) = default;
    SquareMatrix<T>& operator=(SquareMatrix<T>&& other) noexcept = default;

    // Specific Methods for Square Matrixes

    // Convenient getter (since rows == cols)
    int get_size() const;

    // Calculating the Trace of a Matrix (Sum of the Diagonal)
    T trace() const;

    // Transposition (returns a new square matrix)
    SquareMatrix<T>* transpose() const;

    // Symmetry check (A[i][j] == A[j][i])
    [[nodiscard]] bool is_symmetric(double tol = 1e-9) const;

    // Covariant Return for an Interface
    SquareMatrix<T>* clone() const override;

    SquareMatrix<T>* add(const IMatrix<T>& other) const override;
    SquareMatrix<T>* mult(const T& value) const override;
    SquareMatrix<T>* mult(const IMatrix<T>& other) const override;

    // for visibility of all methods from base class Matrix
    using Matrix<T>::operator*;

    // in-place operations
    SquareMatrix<T>& operator+=(const IMatrix<T>& other);
    SquareMatrix<T>& operator-=(const IMatrix<T>& other);
    SquareMatrix<T>& operator*=(const T& scalar);

    // math operations
    SquareMatrix<T> operator+(const IMatrix<T>& other) const;
    SquareMatrix<T> operator-(const IMatrix<T>& other) const;
    SquareMatrix<T> operator*(const T& scalar) const;
    SquareMatrix<T> operator*(const IMatrix<T>& other) const;
};

template <typename T>
SquareMatrix<T> operator*(const T& scalar, const SquareMatrix<T>& m);

#include "square_matrix.tpp"