#include <gtest/gtest.h>
#include "../src/algebra/square_matrix.hpp"
#include "../src/types/complex.hpp"
#include <stdexcept>

// Constructor Testing
TEST(SquareMatrixTest, ConstructorsAndSize) {
    SquareMatrix<int> empty;
    EXPECT_EQ(empty.get_size(), 0);
    EXPECT_EQ(empty.get_rows(), 0);
    EXPECT_EQ(empty.get_cols(), 0);

    SquareMatrix<int> m1(3);
    EXPECT_EQ(m1.get_size(), 3);
    EXPECT_EQ(m1.get_rows(), 3);
    EXPECT_EQ(m1.get_cols(), 3);
    EXPECT_EQ(m1.get(2, 2), 0);

    SquareMatrix<double> m2(2, 3.14);
    EXPECT_DOUBLE_EQ(m2.get(0, 0), 3.14);
    EXPECT_DOUBLE_EQ(m2.get(1, 1), 3.14);

    int data[] = {1, 2, 3, 4}; // 2x2
    SquareMatrix<int> m3(data, 2);
    EXPECT_EQ(m3.get(0, 1), 2);
    EXPECT_EQ(m3.get(1, 0), 3);
}

// Testing Mathematical Overrides (Covariant Return)
TEST(SquareMatrixTest, MathOverrides) {
    int d1[] = {1, 2,
                3, 4};
    int d2[] = {5, 6,
                7, 8};
    SquareMatrix<int> sq1(d1, 2);
    SquareMatrix<int> sq2(d2, 2);

    // 1. Addition
    // The compiler must understand that add returns specifically a SquareMatrix<int>
    SquareMatrix<int>* sum = sq1.add(sq2);
    EXPECT_EQ(sum->get_size(), 2);
    EXPECT_EQ(sum->get(0, 0), 6);
    EXPECT_EQ(sum->get(1, 1), 12);
    delete sum;

    SquareMatrix<int> wrong_sq(3);
    EXPECT_THROW(sq1.add(wrong_sq), std::invalid_argument);

    // 2. Scalar Multiplication
    SquareMatrix<int>* scaled = sq1.mult_scalar(10);
    EXPECT_EQ(scaled->get_size(), 2);
    EXPECT_EQ(scaled->get(0, 1), 20);
    delete scaled;

    // 3. Matrix Multiplication
    SquareMatrix<int>* prod = sq1.mult(sq2);
    EXPECT_EQ(prod->get_size(), 2);
    // [1 2] * [5 6] = [19 22]
    // [3 4]   [7 8]   [43 50]
    EXPECT_EQ(prod->get(0, 0), 19);
    EXPECT_EQ(prod->get(1, 1), 50);
    delete prod;
}

// Testing Specific Methods: Trace
TEST(SquareMatrixTest, Trace) {
    int data[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    SquareMatrix<int> m(data, 3);
    // След: 1 + 5 + 9 = 15
    EXPECT_EQ(m.trace(), 15);

    SquareMatrix<int> empty;
    EXPECT_EQ(empty.trace(), 0);
}

// Testing Specific Methods: Transpose
TEST(SquareMatrixTest, Transpose) {
    int data[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    SquareMatrix<int> m(data, 3);

    SquareMatrix<int>* trans = m.transpose();
    EXPECT_EQ(trans->get_size(), 3);

    EXPECT_EQ(trans->get(0, 1), 4);
    EXPECT_EQ(trans->get(1, 0), 2);
    EXPECT_EQ(trans->get(2, 0), 3);
    EXPECT_EQ(trans->get(2, 2), 9);

    delete trans;
}

    // Testing Specific Methods: is_symmetric
TEST(SquareMatrixTest, IsSymmetric) {
    // Symmetric matrix
    int sym_data[] = {1, 7, 3,
                      7, 4, -5,
                      3, -5, 6};
    SquareMatrix<int> sym_m(sym_data, 3);
    EXPECT_TRUE(sym_m.is_symmetric());

    // Asymmetric matrix
    int asym_data[] = {1, 7, 3,
                       2, 4, -5,  // 2 != 7
                       3, -5, 6};
    SquareMatrix<int> asym_m(asym_data, 3);
    EXPECT_FALSE(asym_m.is_symmetric());

    // Test with Complex numbers and tolerance
    SquareMatrix<Complex<double>> comp_m(2);
    comp_m.set(0, 1, Complex<double>(2.0, 3.0));
    comp_m.set(1, 0, Complex<double>(2.000000001, 3.0));

    // With the default tol = 1e-9, it should return false
    EXPECT_FALSE(comp_m.is_symmetric());
    // With a coarse tol = 1e-6, it should return true
    EXPECT_TRUE(comp_m.is_symmetric(1e-6));
}

// Copy and Move Testing
TEST(SquareMatrixTest, CopyAndMoveSemantics) {
    SquareMatrix<int> original(2, 42);

    SquareMatrix<int> copy = original;
    EXPECT_EQ(copy.get_size(), 2);
    copy.set(0, 0, 99);
    EXPECT_EQ(original.get(0, 0), 42);
    EXPECT_EQ(copy.get(0, 0), 99);

    SquareMatrix<int> moved = std::move(original);
    EXPECT_EQ(moved.get_size(), 2);
    EXPECT_EQ(moved.get(1, 1), 42);
    EXPECT_EQ(original.get_size(), 0);
}