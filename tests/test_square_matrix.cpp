#include <gtest/gtest.h>
#include "../src/algebra/square_matrix.hpp"
#include "../src/types/complex.hpp"
#include <stdexcept>

// Constructor Testing
TEST(SquareMatrixTest, ConstructorsAndSize) {
    SquareMatrix<int> empty;
    EXPECT_EQ(empty.get_size(), 0);

    SquareMatrix<int> m1(3);
    EXPECT_EQ(m1.get_size(), 3);
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
    SquareMatrix<int>* scaled = sq1.mult(10);
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

// Overloaded Operators for SquareMatrix
TEST(SquareMatrixTest, OperatorOverloading) {
    int d1[] = {1, 2, 
                3, 4};
    int d2[] = {5, 6, 
                7, 8};
    SquareMatrix<int> sq1(d1, 2);
    SquareMatrix<int> sq2(d2, 2);

    // 1. In-place operators
    sq1 += sq2;
    EXPECT_EQ(sq1.get(0, 0), 6);
    EXPECT_EQ(sq1.get(1, 1), 12);
    sq1 -= sq2; 

    // 2. Binary Operators (Type checking via compilation)
    SquareMatrix<int> sum = sq1 + sq2;
    EXPECT_EQ(sum.get(0, 0), 6);
    
    SquareMatrix<int> diff = sq1 - sq2;
    EXPECT_EQ(diff.get(0, 1), -4);

    SquareMatrix<int> scaled = sq1 * 10;
    SquareMatrix<int> scaled_global = 10 * sq1;
    EXPECT_EQ(scaled.get(1, 0), 30);
    EXPECT_EQ(scaled_global.get(1, 0), 30);

    // 3. Matrix Multiplication
    SquareMatrix<int> prod = sq1 * sq2;
    EXPECT_EQ(prod.get(0, 0), 19);
    EXPECT_EQ(prod.get(1, 1), 50);

    // 4. Exception check: Multiplying by a rectangular matrix
    Matrix<int> rect(2, 3); // 2x3 matrix
    EXPECT_THROW(sq1 * rect, std::invalid_argument); // could be result 2x3, but SquareMatrix
}

// Inherited Matrix-Vector Multiplication
TEST(SquareMatrixTest, VectorMultiplication) {
    int d1[] = {1, 2, 
                3, 4};
    SquareMatrix<int> sq(d1, 2);
    
    int v_data[] = {1, 1};
    Vector<int> vec(v_data, 2);

    // We didn't write operator* for vector in SquareMatrix, 
    // it inherits from Matrix
    Vector<int> res = sq * vec;
    
    EXPECT_EQ(res.get_size(), 2);
    EXPECT_EQ(res[0], 3); // 1*1 + 2*1
    EXPECT_EQ(res[1], 7); // 3*1 + 4*1
}