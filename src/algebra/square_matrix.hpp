#pragma once

#include "matrix.hpp"
#include "../types/math.hpp"

template <typename T>
class SquareMatrix : public Matrix<T> {
private: 
    //TODO: боремся с прямоугольностью для квадртаных матриц) 
    [[nodiscard]] int get_rows() const override { return Matrix::get_rows(); }
    [[nodiscard]] int get_cols() const override { return Matrix::get_cols(); }
    
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
    SquareMatrix<T>* mult_scalar(const T& value) const override;
    SquareMatrix<T>* mult(const IMatrix<T>& other) const override;
};

#include "square_matrix.tpp"