#include <gtest/gtest.h>
#include "../src/algebra/system_of_equations.hpp"
#include <stdexcept>

// 2x + y = 5
// -x + 3y = 1
// Answer: x = 2, y = 1
SystemOfEquations<double> create_simple_system() {
    double m_data[] = {2.0, 1.0,
                       -1.0, 3.0};
    SquareMatrix<double> A(m_data, 2);

    double b_data[] = {5.0, 1.0};
    Vector<double> b(b_data, 2);

    return SystemOfEquations<double>(A, b);
}

// Testing Constructors and Setters
TEST(SystemOfEquationsTest, ConstructorsAndDimensions) {
    SquareMatrix<double> A(2);
    Vector<double> b(3);

    EXPECT_THROW(SystemOfEquations<double>(A, b), std::invalid_argument);

    SystemOfEquations<double> sys(A, Vector<double>(2));
    EXPECT_EQ(sys.get_size(), 2);

    EXPECT_THROW(sys.set_rhs(Vector<double>(5)), std::invalid_argument);
    EXPECT_THROW(sys.set_matrix(SquareMatrix<double>(5)), std::invalid_argument);
}

// Naive Gaussian Testing
TEST(SystemOfEquationsTest, NaiveGauss) {
    auto sys = create_simple_system();
    Vector<double> x = sys.solve_gauss();

    EXPECT_DOUBLE_EQ(x.get(0), 2.0);
    EXPECT_DOUBLE_EQ(x.get(1), 1.0);

    // Residual check
    EXPECT_NEAR(sys.residual(x), 0.0, 1e-12);
}

    // Gaussian Elimination with Pivoting
TEST(SystemOfEquationsTest, GaussWithPivot_ZeroOnDiagonal) {
    // System:
    // 0x + 1y = 1
    // 2x + 3y = 5
    // Answer: y = 1, x = 1
    double m_data[] = {0.0, 1.0,
                       2.0, 3.0};
    SquareMatrix<double> A(m_data, 2);

    double b_data[] = {1.0, 5.0};
    Vector<double> b(b_data, 2);

    SystemOfEquations<double> sys(A, b);

    // Naive Gauss should fail due to division by zero
    EXPECT_THROW(sys.solve_gauss(), std::runtime_error);

    // Gaussian elimination with pivoting should swap rows and solve successfully
    Vector<double> x = sys.solve_gauss_with_pivot();
    EXPECT_DOUBLE_EQ(x.get(0), 1.0);
    EXPECT_DOUBLE_EQ(x.get(1), 1.0);
}

// PLU Decomposition Testing
TEST(SystemOfEquationsTest, PLU_Decomposition) {
    auto sys = create_simple_system();
    Vector<double> x = sys.solve_plu();

    EXPECT_DOUBLE_EQ(x.get(0), 2.0);
    EXPECT_DOUBLE_EQ(x.get(1), 1.0);
}

// Testing PLU Cache (Change vector b)
TEST(SystemOfEquationsTest, PLU_CacheReuse) {
    auto sys = create_simple_system();

    Vector<double> x1 = sys.solve_plu();
    EXPECT_DOUBLE_EQ(x1.get(0), 2.0);

    // Change b to: {3.0, 2.0}
    // new system:
    // 2x + y = 3
    // -x + 3y = 2
    // Answer: x = 1, y = 1
    double new_b_data[] = {3.0, 2.0};
    sys.set_rhs(Vector<double>(new_b_data, 2));

    Vector<double> x2 = sys.solve_plu();
    EXPECT_DOUBLE_EQ(x2.get(0), 1.0);
    EXPECT_DOUBLE_EQ(x2.get(1), 1.0);
}

// Testing for a Singular (Degenerate) Matrix
TEST(SystemOfEquationsTest, SingularMatrix) {
    // The system has no unique solution (the rows are proportional):
    // 1x + 1y = 2
    // 2x + 2y = 4
    double m_data[] = {1.0, 1.0,
                       2.0, 2.0};
    SquareMatrix<double> A(m_data, 2);

    double b_data[] = {2.0, 4.0};
    Vector<double> b(b_data, 2);

    SystemOfEquations<double> sys(A, b);

    // All methods must detect the singularity and throw an error
    EXPECT_THROW(sys.solve_gauss(), std::runtime_error);
    EXPECT_THROW(sys.solve_gauss_with_pivot(), std::runtime_error);
    EXPECT_THROW(sys.solve_plu(), std::runtime_error);
}

// Testing Metrics and Factories
TEST(SystemOfEquationsTest, MetricsAndFactories) {
    // Generating a 3x3 Hilbert matrix
    auto sys = SystemOfEquations<double>::hilbert(3);
    EXPECT_EQ(sys.get_size(), 3);

    Vector<double> x = sys.solve_plu();

    // The residual must be very small
    EXPECT_NEAR(sys.residual(x), 0.0, 1e-10);

    // Relative Error Test
    Vector<double> exact(3, 1.0);
    Vector<double> approx(3, 0.9);

    // Difference of 0.1 per element
    // Norm of the difference vector: sqrt(0.01 + 0.01 + 0.01) = sqrt(0.03)
    // Norm of the exact vector: sqrt(1 + 1 + 1) = sqrt(3)
    // Error: sqrt(0.03) / sqrt(3) = sqrt(0.01) = 0.1
    EXPECT_NEAR(SystemOfEquations<double>::relative_error(approx, exact), 0.1, 1e-9);
}