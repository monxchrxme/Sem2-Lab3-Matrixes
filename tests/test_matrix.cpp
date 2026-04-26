#include <gtest/gtest.h>
#include "../src/algebra/matrix.hpp"
#include "../src/types/complex.hpp"
#include <stdexcept>

// Testing Constructors and Access
TEST(MatrixTest, ConstructorsAndGetters) {
    Matrix<int> empty;
    EXPECT_EQ(empty.get_rows(), 0);
    EXPECT_EQ(empty.get_cols(), 0);

    Matrix<int> m1(2, 3);
    EXPECT_EQ(m1.get_rows(), 2);
    EXPECT_EQ(m1.get_cols(), 3);
    EXPECT_EQ(m1.get(1, 2), 0);

    Matrix<double> m2(2, 2, 3.14);
    EXPECT_DOUBLE_EQ(m2.get(0, 0), 3.14);
    EXPECT_DOUBLE_EQ(m2.get(1, 1), 3.14);

    int data[] = {1, 2, 3, 4, 5, 6};
    Matrix<int> m3(data, 2, 3);
    EXPECT_EQ(m3.get(0, 0), 1);
    EXPECT_EQ(m3.get(1, 2), 6);
}

// Testing Exceptions for Out-of-Bounds Access
TEST(MatrixTest, OutOfBounds) {
    Matrix<int> m(2, 2);

    EXPECT_THROW(m.get(-1, 0), IndexOutOfRange);
    EXPECT_THROW(m.set(0, -1, 5), IndexOutOfRange);

    EXPECT_THROW(m.get(2, 0), IndexOutOfRange);
    EXPECT_THROW(m.set(0, 2, 5), IndexOutOfRange);
}

// Math Testing (Addition and Scalar Multiplication)
TEST(MatrixTest, AdditionAndScalarMult) {
    int data1[] = {1, 2, 3, 4};
    int data2[] = {5, 6, 7, 8};
    Matrix<int> m1(data1, 2, 2);
    Matrix<int> m2(data2, 2, 2);

    Matrix<int>* sum = m1.add(m2);
    EXPECT_EQ(sum->get(0, 0), 6);
    EXPECT_EQ(sum->get(1, 1), 12);
    delete sum;

    // Adding matrices of different dimensions should throw an exception
    Matrix<int> m_wrong(3, 3);
    EXPECT_THROW(m1.add(m_wrong), std::invalid_argument);

    // Scalar Multiplication
    Matrix<int>* scaled = m1.mult(10);
    EXPECT_EQ(scaled->get(0, 1), 20);
    EXPECT_EQ(scaled->get(1, 0), 30);
    delete scaled;
}

// Matrix Multiplication Testing
TEST(MatrixTest, MatrixMultiplication) {
    // 2x3 * 3x2 = 2x2
    int data1[] = {1, 2, 3,
                   4, 5, 6};
    int data2[] = {7, 8,
                   9, 10,
                   11, 12};
    Matrix<int> m1(data1, 2, 3);
    Matrix<int> m2(data2, 3, 2);

    Matrix<int>* res = m1.mult(m2);

    EXPECT_EQ(res->get_rows(), 2);
    EXPECT_EQ(res->get_cols(), 2);
    // Check element [0][0]: 1*7 + 2*9 + 3*11 = 7 + 18 + 33 = 58
    EXPECT_EQ(res->get(0, 0), 58);
    // Check element [1][1]: 4*8 + 5*10 + 6*12 = 32 + 50 + 72 = 154
    EXPECT_EQ(res->get(1, 1), 154);

    delete res;

    // Multiplication with Mismatched Inner Dimensions
    Matrix<int> m_wrong(4, 4);
    EXPECT_THROW(m1.mult(m_wrong), std::invalid_argument);
}

// Testing Tools for the Gaussian Method
TEST(MatrixTest, GaussRowOperations) {
    int data[] = {1, 2, 3,
                  4, 5, 6,
                  7, 8, 9};
    Matrix<int> m(data, 3, 3);

    // Test swap_rows
    m.swap_rows(0, 2); // Swap the 1st and 3rd rows
    EXPECT_EQ(m.get(0, 0), 7);
    EXPECT_EQ(m.get(2, 0), 1);

    EXPECT_THROW(m.swap_rows(0, 5), IndexOutOfRange);

    // Test row_operation: target = target - factor * source
    // Restore the matrix to its original state
    m.swap_rows(0, 2);
    // We want to zero out element m[1][0] (which is equal to 4)
    // Formula: Row 1 = Row 1 - 4 * Row 0
    m.row_operation(1, 0, 4);

    EXPECT_EQ(m.get(1, 0), 0);  // 4 - 4*1 = 0
    EXPECT_EQ(m.get(1, 1), -3); // 5 - 4*2 = -3
    EXPECT_EQ(m.get(1, 2), -6); // 6 - 4*3 = -6
}

// Frobenius Norm Testing (including Complex)
TEST(MatrixTest, FrobeniusNorm) {
    // simple matrix
    double data[] = {3.0, 4.0, 0.0, 0.0};
    Matrix<double> m1(data, 2, 2);
    // sqrt(3^2 + 4^2 + 0 + 0) = sqrt(25) = 5
    EXPECT_DOUBLE_EQ(m1.norm(), 5.0);

    // matrix with complex numbers
    Matrix<Complex<double>> m2(2, 1);
    m2.set(0, 0, Complex<double>(1.0, 2.0)); // модуль ^ 2 = 5
    m2.set(1, 0, Complex<double>(3.0, 4.0)); // модуль ^ 2 = 25
    // norm = sqrt(5 + 25) = sqrt(30)
    EXPECT_DOUBLE_EQ(m2.norm(), std::sqrt(30.0));
}

// Copy and Move Testing (Rule of Zero)
TEST(MatrixTest, CopyAndMoveSemantics) {
    Matrix<int> original(2, 2, 42);

    // Copy
    Matrix<int> copy = original;
    EXPECT_EQ(copy.get(0, 0), 42);

    // Modifying the copy should not affect the original
    copy.set(0, 0, 99);
    EXPECT_EQ(original.get(0, 0), 42);
    EXPECT_EQ(copy.get(0, 0), 99);

    Matrix<int> moved = std::move(original);
    EXPECT_EQ(moved.get_rows(), 2);
    EXPECT_EQ(moved.get(1, 1), 42);
    // After moving, the original should not point to the old data (DynamicArray behavior)
    EXPECT_EQ(original.get_rows(), 0);
}

// In-place Operators Testing (+=, -=, *=)
TEST(MatrixTest, InPlaceOperators) {
    int d1[] = {1, 2, 3, 4};
    int d2[] = {5, 6, 7, 8};
    Matrix<int> m1(d1, 2, 2);
    Matrix<int> m2(d2, 2, 2);

    // Testing +=
    m1 += m2;
    EXPECT_EQ(m1.get(0, 0), 6);
    EXPECT_EQ(m1.get(1, 1), 12);

    // Testing -=
    m1 -= m2;
    EXPECT_EQ(m1.get(0, 0), 1); 
    EXPECT_EQ(m1.get(1, 1), 4);

    // Testing *= (Scalar)
    m1 *= 10;
    EXPECT_EQ(m1.get(0, 0), 10);
    EXPECT_EQ(m1.get(1, 1), 40);

    // Testing chaining 
    m1 += m2; 
    m1 -= m2; 
    EXPECT_EQ(m1.get(0, 0), 10); 

    // Exception check
    Matrix<int> m_wrong(3, 3);
    EXPECT_THROW(m1 += m_wrong, std::invalid_argument);
    EXPECT_THROW(m1 -= m_wrong, std::invalid_argument);
}

// Binary Arithmetic Operators (+, -, * scalar)
TEST(MatrixTest, BinaryOperators) {
    int d1[] = {10, 20, 30, 40};
    int d2[] = {1, 2, 3, 4};
    Matrix<int> m1(d1, 2, 2);
    Matrix<int> m2(d2, 2, 2);

    // Testing Operator +
    Matrix<int> sum = m1 + m2;
    EXPECT_EQ(sum.get(0, 0), 11);
    EXPECT_EQ(sum.get(1, 1), 44);
    EXPECT_EQ(m1.get(0, 0), 10); 

    // Testing Operator -
    Matrix<int> diff = m1 - m2;
    EXPECT_EQ(diff.get(0, 1), 18);
    EXPECT_EQ(diff.get(1, 0), 27);

    // Testing Scalar Multiplication (Right and Left)
    Matrix<int> scaled_right = m1 * 2;
    Matrix<int> scaled_left = 3 * m1;
    
    EXPECT_EQ(scaled_right.get(0, 0), 20);
    EXPECT_EQ(scaled_left.get(0, 0), 30);
}


// Matrix and Vector Multiplication Operators (*)
TEST(MatrixTest, MultiplicationOperators) {
    // 2x3 Matrix
    int dataA[] = {1, 2, 3,
                   4, 5, 6};
    // 3x2 Matrix
    int dataB[] = {7, 8,
                   9, 10,
                   11, 12};
    Matrix<int> A(dataA, 2, 3);
    Matrix<int> B(dataB, 3, 2);

    // Matrix * Matrix
    Matrix<int> C = A * B;
    EXPECT_EQ(C.get_rows(), 2);
    EXPECT_EQ(C.get_cols(), 2);
    EXPECT_EQ(C.get(0, 0), 58);  // 1*7 + 2*9 + 3*11
    EXPECT_EQ(C.get(1, 1), 154); // 4*8 + 5*10 + 6*12

    // Exception: Inner dimensions mismatch
    Matrix<int> wrong_dim(4, 4);
    EXPECT_THROW(A * wrong_dim, std::invalid_argument);

    // Matrix * Vector
    int data_vec[] = {1, 0, -1};
    Vector<int> v(data_vec, 3);
    
    Vector<int> result_vec = A * v;
    EXPECT_EQ(result_vec.get_size(), 2);
    EXPECT_EQ(result_vec[0], -2); // 1*1 + 2*0 + 3*(-1)
    EXPECT_EQ(result_vec[1], -2); // 4*1 + 5*0 + 6*(-1)

    // Exception: Matrix-Vector dimension mismatch
    Vector<int> wrong_vec(2);
    EXPECT_THROW(A * wrong_vec, std::invalid_argument);
}