#pragma once

#include "imatrix.hpp"
#include "square_matrix.hpp"
#include "../lab2/memory/dynamic_array.hpp"
#include "../lab2/types/exceptions.hpp"
#include <stdexcept>

// Forward Declaration
template <typename T> class Vector;

// Enumeration for Matrix Type
enum class TriangleType {
    Lower, // Lower triangular (zeros above the diagonal)
    Upper  // Upper triangular (zeros below the diagonal)
};

template <typename T>
class TriangularMatrix : public IMatrix<T> {
private:
    DynamicArray<T> m_data; // Store only significant elements
    int m_size;
    TriangleType m_type;

    inline static const T ZERO{};

    bool check_index(int i, int j) const noexcept;
    bool in_structure(int i, int j) const noexcept; // Does it lie in the significant half?
    int flat_index(int i, int j) const noexcept;

public:
    // Constructors
    TriangularMatrix(TriangleType type = TriangleType::Upper);
    explicit TriangularMatrix(int size, TriangleType type = TriangleType::Upper);
    TriangularMatrix(int size, const T& fill_value, TriangleType type = TriangleType::Upper);

    // Rule of Zero (with Custom Move)
    ~TriangularMatrix() override = default;
    TriangularMatrix(const TriangularMatrix<T>& other) = default;
    TriangularMatrix<T>& operator=(const TriangularMatrix<T>& other) = default;
    TriangularMatrix(TriangularMatrix<T>&& other) noexcept;
    TriangularMatrix<T>& operator=(TriangularMatrix<T>&& other) noexcept;

    // getters
    int get_rows() const override;
    int get_cols() const override;
    int get_size() const;
    TriangleType get_type() const;

    // access
    const T& get(int i, int j) const override;
    void set(int i, int j, const T& value) override;

    // math
    SquareMatrix<T>* add(const IMatrix<T>& other) const override;
    TriangularMatrix<T>* mult(const T& value) const override;
    SquareMatrix<T>* mult(const IMatrix<T>& other) const override;

    // In-place math operations
    TriangularMatrix<T>& operator+=(const IMatrix<T>& other);
    TriangularMatrix<T>& operator-=(const IMatrix<T>& other);
    TriangularMatrix<T>& operator*=(const T& scalar);

    // math operations
    SquareMatrix<T> operator+(const IMatrix<T>& other) const;
    SquareMatrix<T> operator-(const IMatrix<T>& other) const;
    TriangularMatrix<T> operator*(const T& scalar) const;
    SquareMatrix<T> operator*(const IMatrix<T>& other) const;
    Vector<T> operator*(const Vector<T>& v) const;

    double norm() const override;
    TriangularMatrix<T>* clone() const override;

    T trace() const;
    T determinant() const;
};

template <typename T>
TriangularMatrix<T> operator*(const T& scalar, const TriangularMatrix<T>& m);

#include "triangular_matrix.tpp"