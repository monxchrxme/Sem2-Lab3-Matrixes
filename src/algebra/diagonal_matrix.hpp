#pragma once

#include "imatrix.hpp"
#include "square_matrix.hpp"
#include "vector.hpp"
#include "../lab2/memory/dynamic_array.hpp"
#include "../lab2/types/exceptions.hpp"
#include "../types/math.hpp"
#include <stdexcept>

// Forward Declaration
template <typename T> class Vector;

template <typename T>
class DiagonalMatrix : public IMatrix<T> {
private:
    DynamicArray<T> m_data;
    int m_size;

    inline static const T ZERO{};

public:
    // Constructors and the Rule of Zero
    DiagonalMatrix();
    explicit DiagonalMatrix(int size);
    DiagonalMatrix(int size, const T& fill_value);
    DiagonalMatrix(T* diag_items, int size);

    ~DiagonalMatrix() override = default;
    DiagonalMatrix(const DiagonalMatrix<T>& other) = default;
    DiagonalMatrix<T>& operator=(const DiagonalMatrix<T>& other) = default;
    DiagonalMatrix(DiagonalMatrix<T>&& other) noexcept;
    DiagonalMatrix<T>& operator=(DiagonalMatrix<T>&& other) noexcept;

    // getters
    int get_rows() const override;
    int get_cols() const override;
    int get_size() const;

    const T& get(int i, int j) const override;
    void set(int i, int j, const T& value) override;

    SquareMatrix<T>* add(const IMatrix<T>& other) const override;
    DiagonalMatrix<T>* mult(const T& value) const override;
    SquareMatrix<T>* mult(const IMatrix<T>& other) const override;

    // In-place math operations
    // if the other has 0's out of diagonal - will be thrown an exception 
    DiagonalMatrix<T>& operator+=(const IMatrix<T>& other);
    DiagonalMatrix<T>& operator-=(const IMatrix<T>& other);
    DiagonalMatrix<T>& operator*=(const T& scalar);

    // math operations (could be 0's out of structure)
    SquareMatrix<T> operator+(const IMatrix<T>& other) const;
    SquareMatrix<T> operator-(const IMatrix<T>& other) const;
    DiagonalMatrix<T> operator*(const T& scalar) const;
    SquareMatrix<T> operator*(const IMatrix<T>& other) const;
    Vector<T> operator*(const Vector<T>& v) const;

    double norm() const override;
    DiagonalMatrix<T>* clone() const override;

    T trace() const;
    T determinant() const;
};

template <typename T>
DiagonalMatrix<T> operator*(const T& scalar, const DiagonalMatrix<T>& m);

#include "diagonal_matrix.tpp"