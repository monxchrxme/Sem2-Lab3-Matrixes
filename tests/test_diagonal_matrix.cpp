#include <gtest/gtest.h>
#include "../src/algebra/diagonal_matrix.hpp"
#include "../src/algebra/square_matrix.hpp"
#include "../src/types/complex.hpp"
#include <stdexcept>
#include <cmath>

// Constructor Testing
TEST(DiagonalMatrixTest, ConstructorsAndSize) {
    // Пустая
    DiagonalMatrix<int> empty;
    EXPECT_EQ(empty.get_size(), 0);

    DiagonalMatrix<int> m1(3);
    EXPECT_EQ(m1.get_size(), 3);
    EXPECT_EQ(m1.get(0, 0), 0);

    DiagonalMatrix<double> m2(4, 2.5);
    EXPECT_DOUBLE_EQ(m2.get(0, 0), 2.5);
    EXPECT_DOUBLE_EQ(m2.get(3, 3), 2.5);

    int data[] = {10, 20, 30};
    DiagonalMatrix<int> m3(data, 3);
    EXPECT_EQ(m3.get(0, 0), 10);
    EXPECT_EQ(m3.get(1, 1), 20);
    EXPECT_EQ(m3.get(2, 2), 30);
}

// Testing Tricky Getters and Setters
TEST(DiagonalMatrixTest, GettersAndSetters) {
    int data[] = {1, 2, 3};
    DiagonalMatrix<int> m(data, 3);

    // Off-diagonal reads must always return 0
    EXPECT_EQ(m.get(0, 1), 0);
    EXPECT_EQ(m.get(2, 0), 0);

    // Writing to the diagonal works
    m.set(1, 1, 99);
    EXPECT_EQ(m.get(1, 1), 99);

    // Writing a non-zero value off the diagonal should throw logic_error
    EXPECT_THROW(m.set(0, 1, 5), std::logic_error);

    // Writing a ZERO off the diagonal should not throw an error (since it does not break the matrix)
    EXPECT_NO_THROW(m.set(0, 1, 0));

    EXPECT_THROW(m.get(3, 3), IndexOutOfRange);
    EXPECT_THROW(m.get(-1, 0), IndexOutOfRange);
}

// Testing Overridden Mathematics
TEST(DiagonalMatrixTest, MathOperations) {
    int d1[] = {1, 2, 3};
    int d2[] = {4, 5, 6};
    DiagonalMatrix<int> diag1(d1, 3);
    DiagonalMatrix<int> diag2(d2, 3);

    // 1. Scalar Multiplication (returns DiagonalMatrix)
    DiagonalMatrix<int>* scaled = diag1.mult(10);
    EXPECT_EQ(scaled->get(0, 0), 10);
    EXPECT_EQ(scaled->get(2, 2), 30);
    EXPECT_EQ(scaled->get(0, 1), 0);
    delete scaled;

    // 2. Addition (returns SquareMatrix)
    SquareMatrix<int>* sum = diag1.add(diag2);
    EXPECT_EQ(sum->get(0, 0), 5); // 1 + 4
    EXPECT_EQ(sum->get(1, 1), 7); // 2 + 5
    EXPECT_EQ(sum->get(0, 1), 0); // 0 + 0
    delete sum;

    // 3. Matrix Multiplication (returns SquareMatrix)
    SquareMatrix<int>* prod = diag1.mult(diag2);
    // When multiplying diagonal matrices, the elements on the main diagonal are multiplied
    EXPECT_EQ(prod->get(0, 0), 4);  // 1 * 4
    EXPECT_EQ(prod->get(1, 1), 10); // 2 * 5
    EXPECT_EQ(prod->get(2, 2), 18); // 3 * 6
    EXPECT_EQ(prod->get(0, 1), 0);
    delete prod;
}

// Dense Matrix Multiplication Test
TEST(DiagonalMatrixTest, MultiplyWithDenseMatrix) {
    int d[] = {2, 3};
    DiagonalMatrix<int> D(d, 2);

    int m_data[] = {1, 2,
                    3, 4};
    SquareMatrix<int> A(m_data, 2);

    SquareMatrix<int>* res = D.mult(A);

    EXPECT_EQ(res->get(0, 0), 2);
    EXPECT_EQ(res->get(0, 1), 4);

    EXPECT_EQ(res->get(1, 0), 9);
    EXPECT_EQ(res->get(1, 1), 12);

    delete res;
}

// Testing Specific O(N) Methods
TEST(DiagonalMatrixTest, SpecificMetrics) {
    int data[] = {2, 3, 4};
    DiagonalMatrix<int> m(data, 3);

    // Trace = 2 + 3 + 4 = 9
    EXPECT_EQ(m.trace(), 9);

    // Determinant = 2 * 3 * 4 = 24
    EXPECT_EQ(m.determinant(), 24);

    // Frobenius's norm = sqrt(2^2 + 3^2 + 4^2) = sqrt(4 + 9 + 16) = sqrt(29)
    EXPECT_DOUBLE_EQ(m.norm(), std::sqrt(29.0));
}

// Copy and Move Testing
TEST(DiagonalMatrixTest, CopyAndMove) {
    int data[] = {5, 5};
    DiagonalMatrix<int> original(data, 2);

    DiagonalMatrix<int> copy = original;
    EXPECT_EQ(copy.get(0, 0), 5);
    copy.set(0, 0, 10);
    EXPECT_EQ(original.get(0, 0), 5);

    DiagonalMatrix<int> moved = std::move(original);
    EXPECT_EQ(moved.get(1, 1), 5);
    EXPECT_EQ(original.get_size(), 0);
}