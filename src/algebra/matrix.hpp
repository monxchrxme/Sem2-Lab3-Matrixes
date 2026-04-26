#pragma once

#include "imatrix.hpp"
#include "../lab2/memory/dynamic_array.hpp"
#include "../lab2/types/exceptions.hpp"

// Forward Declaration
template <typename T> class Vector;

template <typename T>
class Matrix : public IMatrix<T> {
private:
    DynamicArray<T> m_data;
    int m_rows;
    int m_cols;

    // auxiliary private methods (DRY)
    static int checked_size(int rows, int cols);
    //TODO: можно переименовать на check_index для лучшего понимания 
    [[nodiscard]] bool in_bounds(int i, int j) const noexcept; 
    [[nodiscard]] int get_index(int i, int j) const noexcept;

public:
    // constructors
    Matrix(); // create empty matrix 0x0
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, const T& fill_value);
    Matrix(T* items, int rows, int cols);

    // the rule of zero (delegating dynamicArray)
    ~Matrix() override = default;
    Matrix(const Matrix<T>& other) = default;
    Matrix(Matrix<T>&& other) noexcept;
    Matrix<T>& operator=(const Matrix<T>& other) = default;
    Matrix<T>& operator=(Matrix<T>&& other) noexcept;

    [[nodiscard]] int get_rows() const override;
    [[nodiscard]] int get_cols() const override;

    const T& get(int i, int j) const override;
    void set(int i, int j, const T& value) override;

    // covariant return types (returning Matrix* instead of IMatrix*)
    Matrix<T>* add(const IMatrix<T>& other) const override;
    //TODO: можно не писать mult_scalar а оставить просто mult, перегрузка будет
    Matrix<T>* mult_scalar(const T& value) const override;
    Matrix<T>* mult(const IMatrix<T>& other) const override;

    //TODO: добавить операторы именно математические
    // Matrix-vector multiplication (returns a new vector by value)
    Vector<T> operator*(const Vector<T>& v) const;

    [[nodiscard]] double norm() const override;
    Matrix<T>* clone() const override;

    // specific methods for linear algebra
    void swap_rows(int row1, int row2);
    void row_operation(int target_row, int source_row, const T& factor);
    void scale_row(int row, const T& factor);
};

#include "matrix.tpp"