#include <gtest/gtest.h>
#include "../src/algebra/triangular_matrix.hpp"
#include "../src/algebra/square_matrix.hpp"
#include "../src/types/complex.hpp"
#include <stdexcept>
#include <cmath>

// 1. Constructor Testing
TEST(TriangularMatrixTest, ConstructorsAndSize) {
    // Пустая (по умолчанию Upper)
    TriangularMatrix<int> empty;
    EXPECT_EQ(empty.get_size(), 0);
    EXPECT_EQ(empty.get_type(), TriangleType::Upper);

    // With a specified size (Lower)
    TriangularMatrix<int> lower_m(3, TriangleType::Lower);
    EXPECT_EQ(lower_m.get_size(), 3);
    EXPECT_EQ(lower_m.get_type(), TriangleType::Lower);
    EXPECT_EQ(lower_m.get(2, 2), 0);

    // With Fill (Upper)
    TriangularMatrix<double> upper_m(2, 5.5, TriangleType::Upper);
    EXPECT_DOUBLE_EQ(upper_m.get(0, 1), 5.5); // Inside the struct
    EXPECT_DOUBLE_EQ(upper_m.get(1, 0), 0.0); // Outside the structure (must be 0)
}

// 2. Storage and Access Testing (Lower Triangular)
TEST(TriangularMatrixTest, LowerTriangleAccess) {
    TriangularMatrix<int> L(3, TriangleType::Lower);

    // Filling the significant part (below and on the diagonal)
    L.set(0, 0, 1);
    L.set(1, 0, 2);
    L.set(1, 1, 3);
    L.set(2, 0, 4);
    L.set(2, 1, 5);
    L.set(2, 2, 6);

    // Verifying the reading of the significant part
    EXPECT_EQ(L.get(1, 0), 2);
    EXPECT_EQ(L.get(2, 1), 5);

    // Verifying the reading of the zero zone (above the diagonal)
    EXPECT_EQ(L.get(0, 1), 0);
    EXPECT_EQ(L.get(1, 2), 0);

    // Attempt to write a non-zero value to the zero zone
    EXPECT_THROW(L.set(0, 1, 99), std::logic_error);

    // Writing zero to the zero zone (should not fail)
    EXPECT_NO_THROW(L.set(0, 2, 0));
}

// 3. Storage and Access Testing (Upper Triangular)
TEST(TriangularMatrixTest, UpperTriangleAccess) {
    TriangularMatrix<int> U(3, TriangleType::Upper);

    // Fill the significant part (above and on the diagonal)
    U.set(0, 0, 1);
    U.set(0, 1, 2);
    U.set(1, 1, 3);
    U.set(0, 2, 4);
    U.set(1, 2, 5);
    U.set(2, 2, 6);

    // Проверяем чтение значащей части
    EXPECT_EQ(U.get(0, 1), 2);
    EXPECT_EQ(U.get(1, 2), 5);

    // Verify reading of the zero zone (below the diagonal)
    EXPECT_EQ(U.get(1, 0), 0);
    EXPECT_EQ(U.get(2, 0), 0);

    // Attempt to write a non-zero value to the zero zone
    EXPECT_THROW(U.set(2, 0, 99), std::logic_error);
}

// 4. Boundary Testing
TEST(TriangularMatrixTest, OutOfBounds) {
    TriangularMatrix<int> m(2); // 2x2 Upper
    EXPECT_THROW(m.get(-1, 0), IndexOutOfRange);
    EXPECT_THROW(m.set(2, 2, 5), IndexOutOfRange);
}

// 5. Mathematics Testing (Addition and Scalar Multiplication)
TEST(TriangularMatrixTest, AdditionAndScalar) {
    TriangularMatrix<int> U1(2, 2, TriangleType::Upper); // [2 2; 0 2]
    TriangularMatrix<int> U2(2, 3, TriangleType::Upper); // [3 3; 0 3]

    // Scalar Multiplication (returns TriangularMatrix)
    TriangularMatrix<int>* scaled = U1.mult_scalar(10);
    EXPECT_EQ(scaled->get(0, 1), 20);
    EXPECT_EQ(scaled->get(1, 0), 0);
    EXPECT_EQ(scaled->get_type(), TriangleType::Upper);
    delete scaled;

    // Addition (returns SquareMatrix)
    SquareMatrix<int>* sum = U1.add(U2);
    EXPECT_EQ(sum->get(0, 0), 5);
    EXPECT_EQ(sum->get(0, 1), 5);
    EXPECT_EQ(sum->get(1, 0), 0); // 0 + 0 = 0
    delete sum;
}

// 6. Matrix Multiplication Testing
TEST(TriangularMatrixTest, Multiplication) {
    TriangularMatrix<int> L(2, TriangleType::Lower);
    L.set(0, 0, 1);
    L.set(1, 0, 2);
    L.set(1, 1, 3);
    // L = [1 0]
    //     [2 3]

    TriangularMatrix<int> U(2, TriangleType::Upper);
    U.set(0, 0, 4);
    U.set(0, 1, 5);
    U.set(1, 1, 6);
    // U = [4 5]
    //     [0 6]

    // Multiplication L * U (returns SquareMatrix)
    SquareMatrix<int>* prod = L.mult(U);

    // [1*4 + 0*0,  1*5 + 0*6]  = [4, 5]
    // [2*4 + 3*0,  2*5 + 3*6]  = [8, 28]
    EXPECT_EQ(prod->get(0, 0), 4);
    EXPECT_EQ(prod->get(0, 1), 5);
    EXPECT_EQ(prod->get(1, 0), 8);
    EXPECT_EQ(prod->get(1, 1), 28);

    delete prod;
}

// 7. Metrics Testing (O(N) operations)
TEST(TriangularMatrixTest, Metrics) {
    TriangularMatrix<int> U(3, TriangleType::Upper);
    U.set(0, 0, 2);
    U.set(1, 1, 3);
    U.set(2, 2, 4);
    U.set(0, 2, 10); // Off-diagonal element

    // Trace: 2 + 3 + 4 = 9
    EXPECT_EQ(U.trace(), 9);

    // Determinant: 2 * 3 * 4 = 24
    EXPECT_EQ(U.determinant(), 24);

    // Norm: sqrt(2^2 + 3^2 + 4^2 + 10^2) = sqrt(4 + 9 + 16 + 100) = sqrt(129)
    EXPECT_DOUBLE_EQ(U.norm(), std::sqrt(129.0));
}

// 8. Copy and Move Testing
TEST(TriangularMatrixTest, CopyAndMove) {
    TriangularMatrix<int> original(2, 42, TriangleType::Lower);

    TriangularMatrix<int> copy = original;
    EXPECT_EQ(copy.get(1, 0), 42);
    copy.set(1, 0, 99);
    EXPECT_EQ(original.get(1, 0), 42);

    TriangularMatrix<int> moved = std::move(original);
    EXPECT_EQ(moved.get(1, 1), 42);
    EXPECT_EQ(moved.get_type(), TriangleType::Lower);

    EXPECT_EQ(original.get_size(), 0);
}