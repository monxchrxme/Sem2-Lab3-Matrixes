#pragma once

#include "square_matrix.hpp"
#include "triangular_matrix.hpp"
#include "vector.hpp"
#include "../lab2/memory/dynamic_array.hpp"
#include "../types/math.hpp"
#include <stdexcept>
#include <random>
//TODO: include которые не используются в hpp перенести в tpp 

template <typename T>
class SystemOfEquations {
public:
    static constexpr double DEFAULT_TOL = 1e-14;

    // constructors
    SystemOfEquations(const SquareMatrix<T>& A, const Vector<T>& b);

    // getters
    const SquareMatrix<T>& get_matrix() const;
    const Vector<T>& get_rhs() const;
    int get_size() const;

    // setters
    void set(const SquareMatrix<T>& new_A, const Vector<T>& new_b);
    void set_matrix(const SquareMatrix<T>& new_A);
    void set_rhs(const Vector<T>& new_b);

    // SLAU solving methods
    // 1. Naive Gauss (Without Permutations)
    Vector<T> solve_gauss(double tol = DEFAULT_TOL) const;
    // 2. Gauss with partial pivoting
    Vector<T> solve_gauss_with_pivot(double tol = DEFAULT_TOL) const;
    // 3. Naive LU Decomposition (Without Permutations)
    Vector<T> solve_lu(double tol = DEFAULT_TOL) const;
    // 4. Solution via PLU Decomposition (with Caching)
    Vector<T> solve_plu(double tol = DEFAULT_TOL) const;

    // For testing
    static SystemOfEquations<T> random(int n, unsigned seed = 42);
    static SystemOfEquations<T> hilbert(int n);

    // Quality Metrics
    static double relative_error(const Vector<T>& approx, const Vector<T>& exact);
    double residual(const Vector<T>& x) const;

private:
    SquareMatrix<T> A;
    Vector<T> b;

    void decompose_plu(double tol) const;

    mutable bool lu_ready = false;
    mutable DynamicArray<int> P;
    mutable TriangularMatrix<T> L;
    mutable SquareMatrix<T> U;
};

#include "system_of_equations.tpp"