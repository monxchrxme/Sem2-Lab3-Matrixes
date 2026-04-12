#include <gtest/gtest.h>
#include "../src/algebra/vector.hpp"
#include "../src/algebra/matrix.hpp"
#include "../src/types/complex.hpp"
#include <stdexcept>

// Testing Constructors and Access
TEST(VectorTest, ConstructorsAndAccess) {
    Vector<int> empty;
    EXPECT_EQ(empty.get_size(), 0);

    Vector<int> v1(5);
    EXPECT_EQ(v1.get_size(), 5);
    EXPECT_EQ(v1[0], 0);
    EXPECT_EQ(v1[4], 0);

    Vector<double> v2(3, 3.14);
    EXPECT_DOUBLE_EQ(v2[0], 3.14);
    EXPECT_DOUBLE_EQ(v2[2], 3.14);

    int data[] = {10, 20, 30};
    Vector<int> v3(data, 3);
    EXPECT_EQ(v3.get_size(), 3);
    EXPECT_EQ(v3.get(0), 10);
    EXPECT_EQ(v3[1], 20); // test operator[]
    EXPECT_EQ(v3[2], 30);
}

// Boundary and Exception Testing
TEST(VectorTest, OutOfBounds) {
    Vector<int> v(3);

    EXPECT_THROW(v.get(-1), IndexOutOfRange);
    EXPECT_THROW(v.set(-1, 100), IndexOutOfRange);
    EXPECT_THROW(v[-1], IndexOutOfRange);

    EXPECT_THROW(v.get(3), IndexOutOfRange);
    EXPECT_THROW(v.set(5, 100), IndexOutOfRange);
}

// Basic Math Testing (+, -)
TEST(VectorTest, AdditionAndSubtraction) {
    int d1[] = {1, 2, 3};
    int d2[] = {4, 5, 6};
    Vector<int> v1(d1, 3);
    Vector<int> v2(d2, 3);

    Vector<int> sum = v1 + v2;
    EXPECT_EQ(sum[0], 5);
    EXPECT_EQ(sum[1], 7);
    EXPECT_EQ(sum[2], 9);

    Vector<int> diff = v1 - v2;
    EXPECT_EQ(diff[0], -3);
    EXPECT_EQ(diff[1], -3);
    EXPECT_EQ(diff[2], -3);

    // Adding Vectors of Different Lengths
    Vector<int> v_wrong(4);
    EXPECT_THROW(v1 + v_wrong, std::invalid_argument);
}

// Scalar Multiplication Testing
TEST(VectorTest, ScalarMultiplication) {
    int data[] = {1, -2, 3};
    Vector<int> v(data, 3);

    // Vector * Scalar
    Vector<int> res1 = v * 10;
    EXPECT_EQ(res1[0], 10);
    EXPECT_EQ(res1[1], -20);
    EXPECT_EQ(res1[2], 30);

    // Scalar * Vector (Global Operator)
    Vector<int> res2 = 5 * v;
    EXPECT_EQ(res2[0], 5);
    EXPECT_EQ(res2[1], -10);
    EXPECT_EQ(res2[2], 15);
}

// Dot Product Testing
TEST(VectorTest, DotProduct) {
    int d1[] = {1, 2, 3};
    int d2[] = {4, -5, 6};
    Vector<int> v1(d1, 3);
    Vector<int> v2(d2, 3);

    // 1*4 + 2*(-5) + 3*6 = 4 - 10 + 18 = 12
    EXPECT_EQ(v1.dot(v2), 12);

    Vector<int> v_wrong(2);
    EXPECT_THROW(v1.dot(v_wrong), std::invalid_argument);
}

// Testing the Euclidean Norm
TEST(VectorTest, Norm) {
    double data[] = {3.0, 4.0};
    Vector<double> v(data, 2);
    // sqrt(3^2 + 4^2) = sqrt(25) = 5
    EXPECT_DOUBLE_EQ(v.norm(), 5.0);
}

// Matrix-Vector Multiplication Test
TEST(VectorTest, MatrixVectorMultiplication) {
    // Matrix 2x3
    int m_data[] = {1, 2, 3,
                    4, 5, 6};
    Matrix<int> A(m_data, 2, 3);

    // Vector of size 3
    int v_data[] = {1, 0, -1};
    Vector<int> x(v_data, 3);

    // multiplying A * x
    Vector<int> b = A * x;

    // Expected result: A vector of size 2
    // b[0] = 1*1 + 2*0 + 3*(-1) = 1 - 3 = -2
    // b[1] = 4*1 + 5*0 + 6*(-1) = 4 - 6 = -2
    EXPECT_EQ(b.get_size(), 2);
    EXPECT_EQ(b[0], -2);
    EXPECT_EQ(b[1], -2);

    // Exception check for size mismatch
    Vector<int> x_wrong(2);
    EXPECT_THROW(A * x_wrong, std::invalid_argument);
}

// Copy and Move Testing (Rule of Zero)
TEST(VectorTest, CopyAndMoveSemantics) {
    int data[] = {1, 2, 3};
    Vector<int> original(data, 3);

    // Копирование
    Vector<int> copy = original;
    EXPECT_EQ(copy[0], 1);
    copy.set(0, 99);
    EXPECT_EQ(original[0], 1);
    EXPECT_EQ(copy[0], 99);

    // Перемещение
    Vector<int> moved = std::move(original);
    EXPECT_EQ(moved.get_size(), 3);
    EXPECT_EQ(moved[2], 3);
    // We did not implement a custom move operation for the vector; therefore, the get_size() of the original
    // object may depend on the specific implementation of DynamicArray. However, the crucial point is that no data is leaked
    EXPECT_EQ(original.get_size(), 0);
}