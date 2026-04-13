#pragma once
#include "system_of_equations.hpp"

template <typename T>
SystemOfEquations<T>::SystemOfEquations(const SquareMatrix<T>& matrix, const Vector<T>& rhs) :
    A(matrix),
    b(rhs),
    lu_ready(false),
    P(0),
    L(TriangleType::Lower),
    U(0) {
    if (A.get_size() != b.get_size()) {
        throw std::invalid_argument("SystemOfEquations: Sizes must match");
    }
}

template <typename T>
const SquareMatrix<T>& SystemOfEquations<T>::get_matrix() const { return A; }

template <typename T>
const Vector<T>& SystemOfEquations<T>::get_rhs() const { return b; }

template <typename T>
int SystemOfEquations<T>::get_size() const { return A.get_size(); }

template <typename T>
void SystemOfEquations<T>::set(const SquareMatrix<T>& new_A, const Vector<T>& new_b) {
    if (new_A.get_size() != new_b.get_size()) {
        throw std::invalid_argument("Sizes mismatch");
    }
    A = new_A;
    b = new_b;
    lu_ready = false;
}

template <typename T>
void SystemOfEquations<T>::set_matrix(const SquareMatrix<T>& new_A) {
    if (new_A.get_size() != b.get_size()) {
        throw std::invalid_argument("Sizes mismatch");
    }
    A = new_A;
    lu_ready = false;
}

template <typename T>
void SystemOfEquations<T>::set_rhs(const Vector<T>& new_b) {
    if (A.get_size() != new_b.get_size()) {
        throw std::invalid_argument("Sizes mismatch");
    }
    b = new_b;
}

// Naive Gauss
template <typename T>
Vector<T> SystemOfEquations<T>::solve_gauss(double tol) const {
    int n = get_size();
    SquareMatrix<T> tA = A;
    Vector<T> tb = b;

    for (int k = 0; k < n - 1; ++k) {
        if (static_cast<double>(custom_math::abs(tA.get(k, k))) <= tol) {
            throw std::runtime_error("Naive Gauss: Zero pivot");
        }
        for (int i = k + 1; i < n; ++i) {
            T factor = tA.get(i, k) / tA.get(k, k);
            tA.row_operation(i, k, factor);
            tb.set(i, tb.get(i) - factor * tb.get(k));
        }
    }

    if (static_cast<double>(custom_math::abs(tA.get(n - 1, n - 1))) <= tol) {
        throw std::runtime_error("Naive Gauss: Matrix is singular");
    }

    Vector<T> x(n);
    for (int i = n - 1; i >= 0; --i) {
        T sum = T{};
        for (int j = i + 1; j < n; ++j) {
            sum = sum + tA.get(i, j) * x[j];
        }
        x.set(i, (tb.get(i) - sum) / tA.get(i, i));
    }
    return x;
}

// Gauss with pivot
template <typename T>
Vector<T> SystemOfEquations<T>::solve_gauss_with_pivot(double tol) const {
    int n = get_size();

    SquareMatrix<T> tA = A;
    Vector<T> tb = b;

    // Forward Elimination
    for (int k = 0; k < n - 1; ++k) {
        // 1. Finding Pivot
        int pivot_row = k;
        double max_val = static_cast<double>(custom_math::abs(tA.get(k, k)));

        for (int i = k + 1; i < n; ++i) {
            double val = static_cast<double>(custom_math::abs(tA.get(i, k)));
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }

        if (max_val <= tol) {
            throw std::runtime_error("Gauss with pivot: Zero or near-zero pivot");
        }

        // 2. Permutation of rows in matrix tA and vector tb
        if (pivot_row != k) {
            tA.swap_rows(k, pivot_row);

            T temp_b = tb.get(k);
            tb.set(k, tb.get(pivot_row));
            tb.set(pivot_row, temp_b);
        }

        // 3. Zeroing out elements below the main diagonal
        for (int i = k + 1; i < n; ++i) {
            T factor = tA.get(i, k) / tA.get(k, k);
            tA.row_operation(i, k, factor);

            tb.set(i, tb.get(i) - factor * tb.get(k));
        }
    }

    if (static_cast<double>(custom_math::abs(tA.get(n - 1, n - 1))) <= tol) {
        throw std::runtime_error("Gauss with pivot: Matrix is singular");
    }

    // Backward Substitution
    Vector<T> x(n);
    for (int i = n - 1; i >= 0; --i) {
        T sum = T{};
        for (int j = i + 1; j < n; ++j) {
            sum = sum + tA.get(i, j) * x.get(j);
        }
        x.set(i, (tb.get(i) - sum) / tA.get(i, i));
    }

    return x;
}

// PLU decomposition
template <typename T>
void SystemOfEquations<T>::decompose_plu(double tol) const {
    if (lu_ready) return;
    int n = get_size();
    P = DynamicArray<int>(n);
    for (int i = 0; i < n; ++i) P.set(i, i);
    U = A;
    L = TriangularMatrix<T>(n, TriangleType::Lower);
    for (int i = 0; i < n; ++i) {
        L.set(i, i, T(1));
    }

    for (int k = 0; k < n - 1; ++k) {
        int p_row = k;
        double max_v = static_cast<double>(custom_math::abs(U.get(k, k)));
        for (int i = k + 1; i < n; ++i) {
            double v = static_cast<double>(custom_math::abs(U.get(i, k)));
            if (v > max_v) {
                max_v = v; p_row = i;
            }
        }
        if (max_v <= tol) {
            throw std::runtime_error("Matrix is singular");
        }
        if (p_row != k) {
            U.swap_rows(k, p_row);
            int tp = P.get(k); P.set(k, P.get(p_row)); P.set(p_row, tp);
            for (int j = 0; j < k; ++j) {
                T tl = L.get(k, j); L.set(k, j, L.get(p_row, j)); L.set(p_row, j, tl);
            }
        }
        for (int i = k + 1; i < n; ++i) {
            T f = U.get(i, k) / U.get(k, k);
            L.set(i, k, f);
            for (int j = k; j < n; ++j) {
                U.set(i, j, U.get(i, j) - f * U.get(k, j));
            }
        }
    }

    if (static_cast<double>(custom_math::abs(U.get(n - 1, n - 1))) <= tol) {
        throw std::runtime_error("PLU Decomposition: Matrix is singular");
    }

    lu_ready = true;
}

// Solving via PLU
template <typename T>
Vector<T> SystemOfEquations<T>::solve_plu(double tol) const {
    int n = get_size();
    decompose_plu(tol);
    Vector<T> Pb(n);
    for (int i = 0; i < n; ++i) {
        Pb.set(i, b.get(P.get(i)));
    }

    Vector<T> y(n);
    for (int i = 0; i < n; ++i) {
        T s = T{};
        for (int j = 0; j < i; ++j) {
            s = s + L.get(i, j) * y[j];
        }
        y.set(i, (Pb[i] - s) / L.get(i, i));
    }

    Vector<T> x(n);
    for (int i = n - 1; i >= 0; --i) {
        T s = T{};
        for (int j = i + 1; j < n; ++j) {
            s = s + U.get(i, j) * x[j];
        }
        x.set(i, (y[i] - s) / U.get(i, i));
    }
    return x;
}

// Metrics and fabrics
template <typename T>
double SystemOfEquations<T>::residual(const Vector<T>& x) const {
    return (A * x - b).norm();
}

template <typename T>
double SystemOfEquations<T>::relative_error(const Vector<T>& approx, const Vector<T>& exact) {
    double en = exact.norm();
    return (en == 0.0) ? (approx - exact).norm() : (approx - exact).norm() / en;
}

template <typename T>
SystemOfEquations<T> SystemOfEquations<T>::hilbert(int n) {
    SquareMatrix<T> H(n);
    Vector<T> b(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H.set(i, j, T(1) / T(i + j + 1));
        }
        b.set(i, T(1));
    }
    return SystemOfEquations<T>(H, b);
}

template <typename T>
SystemOfEquations<T> SystemOfEquations<T>::random(int n, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    SquareMatrix<T> A(n);
    Vector<T> b(n);
    for (int i = 0; i < n; ++i) {
        b.set(i, T(dist(gen)));
        for (int j = 0; j < n; ++j) {
            A.set(i, j, T(dist(gen)));
        }
    }
    return SystemOfEquations<T>(A, b);
}