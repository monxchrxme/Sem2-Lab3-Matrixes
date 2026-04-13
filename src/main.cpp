#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <random>
#include <stdexcept>
#include <chrono>
#include <type_traits> // Для генерации случайных чисел в зависимости от типа

//#include "lab2/sequences/array_sequence.hpp"
#include "lab2/sequences/mutable_array_sequence.hpp"

#include "algebra/imatrix.hpp"
#include "algebra/matrix.hpp"
#include "algebra/square_matrix.hpp"
#include "algebra/vector.hpp"
#include "algebra/triangular_matrix.hpp"
#include "algebra/diagonal_matrix.hpp"
#include "algebra/system_of_equations.hpp"
#include "algebra/matrix_vector_io.hpp"
#include "types/complex.hpp"

template <typename T>
using Seq = MutableArraySequence<IMatrix<T>*>;

// ===================================================================
// Auxiliary Input Functions (Safe Reading)
// ===================================================================

static void clear_input() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

static int read_int(const std::string& prompt) {
    int val;
    for (;;) {
        std::cout << prompt;

        std::cin >> val;

        if (!std::cin.fail()) {
            return val;
        }

        if (std::cin.eof()) {
            throw std::runtime_error("stdin closed");
        }

        std::cout << " [!] Invalid input, try again.\n";
        clear_input();
    }
}

template <typename T>
static T read_value(const std::string& prompt) {
    T val;
    for (;;) {
        std::cout << prompt;

        std::cin >> val;

        if (!std::cin.fail()) {
            return val;
        }

        if (std::cin.eof()) {
            throw std::runtime_error("stdin closed");
        }

        std::cout << " [!] Invalid input, try again.\n";
        clear_input();
    }
}

template <typename T>
static int read_valid_idx(const Seq<T>& seq) {
    for (;;) {
        int idx = read_int(" Matrix index: ");
        if (idx >= 0 && idx < seq.get_length()) {
            return idx;
        }
        std::cout << " [!] Index out of range [0, " << seq.get_length() - 1 << "].\n";
    }
}

template <typename T>
static T generate_random(std::mt19937& rng, std::uniform_real_distribution<double>& dist) {
    if constexpr (std::is_same_v<T, double>) {
        return dist(rng);
    } else {
        return T(dist(rng), dist(rng));
    }
}

// ===================================================================
// Matrix Utilities
// ===================================================================

template <typename T>
static const char* matrix_type_name(const IMatrix<T>* m) {
    if (dynamic_cast<const DiagonalMatrix<T>*>(m))   return "Diagonal  ";
    if (dynamic_cast<const TriangularMatrix<T>*>(m)) return "Triangular";
    if (dynamic_cast<const SquareMatrix<T>*>(m))     return "Square    ";
    if (dynamic_cast<const Matrix<T>*>(m))           return "Dense     ";
    return "Unknown   ";
}

template <typename T>
static void show_all(const Seq<T>& seq) {
    if (seq.get_length() == 0) {
        std::cout << " (empty)\n";
        return;
    }
    for (int i = 0; i < seq.get_length(); i++) {
        const IMatrix<T>* m = seq.get(i);
        std::cout << " [" << i << "] "
                  << matrix_type_name(m) << "  "
                  << m->get_rows() << "x" << m->get_cols()
                  << "  norm = " << std::fixed << std::setprecision(4) << m->norm()
                  << "\n";
    }
}

template <typename T>
static void fill_manually(IMatrix<T>* m) {
    std::cout << " Enter elements row by row:\n";
    for (int i = 0; i < m->get_rows(); i++) {
        for (int j = 0; j < m->get_cols(); j++) {
            T v = read_value<T>(" [" + std::to_string(i) + "][" + std::to_string(j) + "] = ");
            try {
                m->set(i, j, v);
            }
            catch (const std::invalid_argument&) {
                std::cout << " [!] Position (" << i << "," << j << ") outside structure, skipped.\n";
            }
            catch (const std::logic_error&) {
                std::cout << " [!] Position (" << i << "," << j << ") is restricted by matrix type, skipped.\n";
            }
        }
    }
}

template <typename T>
static void fill_randomly(IMatrix<T>* m, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < m->get_rows(); i++) {
        for (int j = 0; j < m->get_cols(); j++) {
            try {
                m->set(i, j, generate_random<T>(rng, dist));
            }
            catch (...) {} // Ignore write attempts outside the structure of diagonal/triangular matrices
        }
    }
}

template <typename T>
static IMatrix<T>* create_matrix() {
    std::cout << "\n Matrix type:\n"
              << "  1. Matrix (rectangular dense)\n"
              << "  2. SquareMatrix\n"
              << "  3. DiagonalMatrix\n"
              << "  4. TriangularMatrix Lower\n"
              << "  5. TriangularMatrix Upper\n";
    int type = read_int(" Choice: ");

    IMatrix<T>* m = nullptr;

    try {
        switch (type) {
            case 1: m = new Matrix<T>(read_int(" Rows: "), read_int(" Cols: ")); break;
            case 2: m = new SquareMatrix<T>(read_int(" Size n: "));              break;
            case 3: m = new DiagonalMatrix<T>(read_int(" Size n: "));            break;
            case 4: m = new TriangularMatrix<T>(read_int(" Size n: "), TriangleType::Lower); break;
            case 5: m = new TriangularMatrix<T>(read_int(" Size n: "), TriangleType::Upper); break;
            default: std::cout << " [!] Unknown type.\n"; return nullptr;
        }
    } catch (const std::invalid_argument& e) {
        std::cout << " [!] " << e.what() << "\n";
        return nullptr;
    }

    std::cout << " Fill: 1. Manual  2. Random\n";
    int fill_choice = read_int(" Choice: ");

    if (fill_choice == 1) {
        fill_manually(m);
    } else {
        fill_randomly(m, static_cast<unsigned>(read_int(" Seed: ")));
    }

    return m;
}

template <typename T>
static void matrix_ops(Seq<T>& seq, int idx) {
    for (;;) {
        IMatrix<T>* m = seq.get(idx);
        std::cout << "\n--- Matrix [" << idx << "] "
                  << matrix_type_name(m) << " "
                  << m->get_rows() << "x" << m->get_cols() << " ---\n"
                  << " 1. Print\n"
                  << " 2. get(i, j)\n"
                  << " 3. set(i, j, val)\n"
                  << " 4. norm()\n"
                  << " 5. add(other) -> save result\n"
                  << " 6. mult_scalar(val) -> save result\n"
                  << " 7. mult(other) -> save result\n"
                  << " 0. Back\n";

        int op = read_int(" Choice: ");
        if (op == 0) return;

        switch (op) {
            case 1:
                std::cout << *m;
                break;
            case 2: {
                int i = read_int(" Row i: ");
                int j = read_int(" Col j: ");
                try {
                    std::cout << " [" << i << "][" << j << "] = " << m->get(i, j) << "\n";
                }
                catch (const std::exception& e) {
                    std::cout << " [!] " << e.what() << "\n";
                }
                break;
            }
            case 3: {
                int i = read_int(" Row i: ");
                int j = read_int(" Col j: ");
                T v = read_value<T>(" Value: ");
                try {
                    m->set(i, j, v);
                    std::cout << " Done.\n";
                }
                catch (const std::exception& e) {
                    std::cout << " [!] " << e.what() << "\n";
                }
                break;
            }
            case 4:
                std::cout << " norm = " << std::fixed << std::setprecision(8) << m->norm() << "\n";
                break;
            case 5: {
                show_all(seq);
                int other = read_valid_idx(seq);
                try {
                    IMatrix<T>* res = m->add(*seq.get(other));
                    seq.append(res);
                    std::cout << " Result saved.\n";
                } catch (const std::exception& e) {
                    std::cout << " [!] " << e.what() << "\n";
                }
                break;
            }
            case 6: {
                T sc = read_value<T>(" Scalar: ");
                try {
                    IMatrix<T>* res = m->mult_scalar(sc);
                    seq.append(res);
                    std::cout << " Result saved.\n";
                } catch (const std::exception& e) {
                    std::cout << " [!] " << e.what() << "\n";
                }
                break;
            }
            case 7: {
                show_all(seq);
                int other = read_valid_idx(seq);
                try {
                    IMatrix<T>* res = m->mult(*seq.get(other));
                    seq.append(res);
                    std::cout << " Result saved.\n";
                } catch (const std::exception& e) {
                    std::cout << " [!] " << e.what() << "\n";
                }
                break;
            }
            default:
                std::cout << " [!] Unknown command.\n";
        }
    }
}

template <typename T>
static void run_matrix_menu(Seq<T>& seq) {
    for (;;) {
        std::cout << "\n=== Matrixes (" << seq.get_length() << " stored) ===\n"
                  << " 1. Show all\n"
                  << " 2. Create & save\n"
                  << " 3. Select & operate\n"
                  << " 4. Delete\n"
                  << " 0. Back\n";

        int choice = read_int("Choice: ");
        if (choice == 0) return;

        switch (choice) {
            case 1:
                show_all(seq);
                break;
            case 2: {
                IMatrix<T>* m = create_matrix<T>();
                if (m) {
                    seq.append(m);
                    std::cout << " Saved.\n";
                }
                break;
            }
            case 3: {
                if (seq.get_length() == 0) {
                    std::cout << " No matrixes.\n";
                    break;
                }
                show_all(seq);
                matrix_ops(seq, read_valid_idx(seq));
                break;
            }
            case 4: {
                if (seq.get_length() == 0) {
                    std::cout << " No matrixes.\n";
                    break;
                }
                show_all(seq);
                int idx = read_valid_idx(seq);
                delete seq.get(idx);

                seq.remove_at(idx);

                std::cout << " Deleted. Remaining: " << seq.get_length() << "\n";
                break;
            }
            default:
                std::cout << " [!] Unknown command.\n";
        }
    }
}

// ===================================================================
// Solving Systems of Linear Equations and Experiments
// ===================================================================

template <typename T>
static void print_solution(const Vector<T>& x, int n) {
    for (int i = 0; i < n; i++) {
        std::cout << " x[" << i << "] = "
                  << std::fixed << std::setprecision(8) << x.get(i) << "\n";
    }
}

template <typename T>
static void run_soe_work() {
    int n = read_int("\n Number of equations n: ");
    if (n <= 0) {
        std::cout << " [!] n must be positive.\n";
        return;
    }

    if (n > 1000) {
        std::cout << " [!] Warning: n > 1000 may take a long time.\n";
    }

    SystemOfEquations<T>* soe = nullptr;
    std::cout << " Fill: 1. Manual  2. Random\n";
    int choice = read_int(" Choice: ");

    if (choice == 2) {
        unsigned seed = static_cast<unsigned>(read_int(" Seed: "));
        SystemOfEquations<T> tmp = SystemOfEquations<T>::random(n, seed);
        soe = new SystemOfEquations<T>(tmp);
    } else {
        SquareMatrix<T> A(n);
        Vector<T> b(n);
        std::cout << " Enter matrix A: (if Complex Type use format (a, b))\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                T val = read_value<T>(" A[" + std::to_string(i) + "][" + std::to_string(j) + "] = ");
                A.set(i, j, val);
            }
        }
        std::cout << " Enter vector b:\n";
        for (int i = 0; i < n; i++) {
            T val = read_value<T>(" b[" + std::to_string(i) + "] = ");
            b.set(i, val);
        }
        soe = new SystemOfEquations<T>(A, b);
    }

    std::cout << "\nMatrix A:\n" << soe->get_matrix();
    std::cout << "Vector b:\n" << soe->get_rhs() << "\n";

    for (;;) {
        std::cout << "\n--- SOE (" << n << "x" << n << ") ---\n"
                  << " 1. Solve - Naive Gauss\n"
                  << " 2. Solve - Gauss + pivot\n"
                  << " 3. Solve - Pure LU (No pivot)\n"
                  << " 4. Solve - PLU (Advanced)\n"
                  << " 5. Print system\n"
                  << " 0. Back\n";

        int op = read_int(" Choice: ");
        if (op == 0) break;
        if (op == 5) {
            std::cout << "Matrix A:\n" << soe->get_matrix();
            std::cout << "Vector b:\n" << soe->get_rhs() << "\n";
            continue;
        }
        if (op < 1 || op > 4) {
            std::cout << " [!] Unknown command.\n";
            continue;
        }

        try {
            Vector<T> x;
            if      (op == 1) x = soe->solve_gauss();
            else if (op == 2) x = soe->solve_gauss_with_pivot();
            else if (op == 3) x = soe->solve_lu();
            else              x = soe->solve_plu();

            std::cout << " Solution:\n";
            print_solution(x, n);
            std::cout << " Residual ||Ax-b|| = "
                      << std::scientific << std::setprecision(4)
                      << soe->residual(x) << "\n";
        } catch (const std::exception& e) {
            std::cout << " [!] " << e.what() << "\n";
        }
    }
    delete soe;
}

using Clock = std::chrono::high_resolution_clock;
using Ms    = std::chrono::duration<double, std::milli>;

template <typename T>
static void experiment_4_1() {
    const int sizes[] = {100, 200, 500};
    const unsigned seed = 67;

    std::cout << "\n--- Experiment 4.1: Time for one system (seed=" << seed << ") ---\n";
    std::cout << std::setw(6)  << "n"
              << std::setw(12) << "Gauss"
              << std::setw(12) << "G+pivot"
              << std::setw(12) << "Pure LU"
              << std::setw(16) << "PLU decomp"
              << std::setw(14) << "PLU subst"
              << "\n" << std::string(72, '-') << "\n";

    for (int n : sizes) {
        SystemOfEquations<T> soe = SystemOfEquations<T>::random(n, seed);

        auto t0 = Clock::now();
        Vector<T> x_gauss = soe.solve_gauss();
        double t_gauss = Ms(Clock::now() - t0).count();

        soe = SystemOfEquations<T>::random(n, seed);
        t0  = Clock::now();
        Vector<T> x_piv = soe.solve_gauss_with_pivot();
        double t_pivot = Ms(Clock::now() - t0).count();

        soe = SystemOfEquations<T>::random(n, seed);
        t0  = Clock::now();
        Vector<T> x_lu = soe.solve_lu();
        double t_lu = Ms(Clock::now() - t0).count();

        SystemOfEquations<T> soe_plu = SystemOfEquations<T>::random(n, seed);
        t0 = Clock::now();
        Vector<T> x_full = soe_plu.solve_plu();
        double t_full = Ms(Clock::now() - t0).count();

        t0 = Clock::now();
        Vector<T> x_subst = soe_plu.solve_plu();
        double t_subst = Ms(Clock::now() - t0).count();

        std::cout << std::setw(6)  << n
                  << std::setw(12) << std::fixed << std::setprecision(2) << t_gauss
                  << std::setw(12) << t_pivot
                  << std::setw(12) << t_lu
                  << std::setw(16) << (t_full - t_subst)
                  << std::setw(14) << t_subst
                  << "\n";
    }
}

template <typename T>
static void experiment_4_2() {
    const int n = 250;
    const int ks[] = {1, 10, 50};
    const unsigned seed = 67;

    std::cout << "\n--- Experiment 4.2: Multiple RHS (n=" << n << ", seed=" << seed << ") ---\n";
    std::cout << std::setw(6)  << "k"
              << std::setw(24) << "Gauss+piv k*x (ms)"
              << std::setw(24) << "PLU decomp + k*sub (ms)"
              << "\n" << std::string(54, '-') << "\n";

    std::mt19937 rng(seed + 1);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    MutableArraySequence<Vector<T>> rhs;
    for (int i = 0; i < 50; i++) {
        Vector<T> v(n);
        for (int j = 0; j < n; j++) {
            v.set(j, generate_random<T>(rng, dist));
        }
        rhs.append(v);
    }

    for (int k : ks) {
        SystemOfEquations<T> soe_g = SystemOfEquations<T>::random(n, seed);
        auto t0 = Clock::now();
        for (int i = 0; i < k; i++) {
            soe_g.set_rhs(rhs.get(i));
            Vector<T> x = soe_g.solve_gauss_with_pivot();
        }
        double t_gauss = Ms(Clock::now() - t0).count();

        SystemOfEquations<T> soe_plu = SystemOfEquations<T>::random(n, seed);
        t0 = Clock::now();
        for (int i = 0; i < k; i++) {
            soe_plu.set_rhs(rhs.get(i));
            Vector<T> x = soe_plu.solve_plu();
        }
        double t_plu = Ms(Clock::now() - t0).count();

        std::cout << std::setw(6)  << k
                  << std::setw(24) << std::fixed << std::setprecision(2) << t_gauss
                  << std::setw(24) << t_plu
                  << "\n";
    }
}

template <typename T>
static void experiment_4_3(double tol) {
    const int sizes[] = {5, 10, 15};

    std::cout << "\n--- Experiment 4.3: Hilbert accuracy (x_exact = [1,...,1]) ---\n";

    for (int n : sizes) {
        std::cout << "\n n = " << n << ":\n";
        std::cout << " " << std::setw(14) << "Method"
                  << std::setw(16) << "Residual"
                  << std::setw(16) << "Rel.Error"
                  << "\n " << std::string(44, '-') << "\n";

        SystemOfEquations<T> soe = SystemOfEquations<T>::hilbert(n);
        Vector<T> exact(n);
        for (int i = 0; i < n; i++) {
            exact.set(i, T(1.0));
        }

        const char* names[] = {"Gauss", "Gauss+pivot", "Pure LU", "PLU"};
        for (int m = 0; m < 4; m++) {
            try {
                Vector<T> x;
                if      (m == 0) x = soe.solve_gauss();
                else if (m == 1) x = soe.solve_gauss_with_pivot(tol);
                else if (m == 2) x = soe.solve_lu(tol);
                else             x = soe.solve_plu(tol);

                std::cout << " " << std::setw(14) << names[m]
                          << std::setw(16) << std::scientific << std::setprecision(4)
                          << soe.residual(x)
                          << std::setw(16) << SystemOfEquations<T>::relative_error(x, exact)
                          << "\n";
            } catch (const std::exception& e) {
                std::cout << " " << std::setw(14) << names[m]
                          << " [!] " << e.what() << "\n";
            }
        }
    }
}

template <typename T>
static void run_experiments() {
    std::cout << "\n=== Running all experiments ===\n";
    try {
        experiment_4_1<T>();
    }
    catch (const std::exception& e) {
        std::cout << "[!] Exp 4.1: " << e.what() << "\n";
    }
    try {
        experiment_4_2<T>();
    }
    catch (const std::exception& e) {
        std::cout << "[!] Exp 4.2: " << e.what() << "\n";
    }
    try {
        experiment_4_3<T>(1e-14);
    }
    catch (const std::exception& e) {
        std::cout << "[!] Exp 4.3: " << e.what() << "\n";
    }
    std::cout << "=== Done ===\n";
}

template <typename T>
static void run_soe_menu() {
    for (;;) {
        std::cout << "\n=== System of Linear Equations ===\n"
                  << " 1. Work\n"
                  << " 2. Experiments\n"
                  << " 0. Back\n";
        int c = read_int("Choice: ");
        if (c == 0) return;

        if      (c == 1) run_soe_work<T>();
        else if (c == 2) run_experiments<T>();
        else std::cout << " [!] Unknown command.\n";
    }
}

template <typename T>
static void run_app() {
    Seq<T> matrixes;

    for (;;) {
        std::cout << "\n==============================\n"
                  << "         MAIN MENU\n"
                  << "==============================\n"
                  << " 1. Matrixes\n"
                  << " 2. System of Linear Equations\n"
                  << " 0. Back to Data Type Selection\n";

        int c = read_int("Choice: ");
        if      (c == 0) break;
        else if (c == 1) run_matrix_menu<T>(matrixes);
        else if (c == 2) run_soe_menu<T>();
        else std::cout << " [!] Unknown command.\n";
    }

    for (int i = 0; i < matrixes.get_length(); i++) {
        delete matrixes.get(i);
    }
}

int main() {
    for (;;) {
        std::cout << "\n==============================\n"
                  << "      CHOOSE DATA TYPE\n"
                  << "==============================\n"
                  << " 1. Real numbers (double)\n"
                  << " 2. Complex numbers (Complex)\n"
                  << " 0. Exit\n";

        int c = read_int("Choice: ");
        if      (c == 0) break;
        else if (c == 1) run_app<double>();
        else if (c == 2) run_app<Complex<double>>();
        else std::cout << " [!] Unknown command.\n";
    }

    return 0;
}