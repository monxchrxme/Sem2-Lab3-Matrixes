#pragma once 

#include "imatrix.hpp"
#include "../lab2/memory/dynamic_array.hpp"
#include "../types/math.hpp" 
#include <stdexcept>
#include <cmath>

template <typename T>
class SparseMatrix : public IMatrix<T> {
private:
    struct Triplet {
        int row;
        int col;
        T val;
    };

    int m_rows;
    int m_cols;
    DynamicArray<Triplet> m_data;

    inline static const T ZERO{};

    int find_index(int i, int j) const;

public: 
    SparseMatrix(int rows, int cols);

    SparseMatrix(const SparseMatrix<T>& other) = default;
    SparseMatrix(SparseMatrix<T>&& other) noexcept = default;
    SparseMatrix<T>& operator=(const SparseMatrix<T>& other) = default;
    SparseMatrix<T>& operator=(SparseMatrix<T>&& other) noexcept = default;
    ~SparseMatrix() override = default; 

    [[nodiscard]] int get_rows() const override;
    [[nodiscard]] int get_cols() const override;

    const T& get(int i, int j) const override;
    void set(int i, int j, const T& value) override;

    [[nodiscard]] double norm() const override;
    SparseMatrix<T>* clone() const override;

    SparseMatrix<T>* add(const IMatrix<T>& other) const override;
    SparseMatrix<T>* mult(const T& value) const override;
    SparseMatrix<T>* mult(const IMatrix<T>& other) const override;
    
    bool check_index(int i, int j) const noexcept;

    int get_non_zero_count() const {
        return m_data.get_size();
    }
};

#include "sparse_matrix.tpp"