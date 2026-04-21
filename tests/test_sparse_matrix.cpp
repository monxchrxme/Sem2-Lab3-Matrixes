#include <gtest/gtest.h>
#include "../src/algebra/sparse_matrix.hpp"
#include "../src/types/complex.hpp" 

TEST(SparseMatrixTest, InitializationAndDimensions) {
    SparseMatrix<double> m(5, 10);
    EXPECT_EQ(m.get_rows(), 5);
    EXPECT_EQ(m.get_cols(), 10);
    EXPECT_EQ(m.get_non_zero_count(), 0); 
    
    EXPECT_THROW(SparseMatrix<double>(0, 5), std::invalid_argument);
    EXPECT_THROW(SparseMatrix<double>(5, -2), std::invalid_argument);
}

TEST(SparseMatrixTest, SetAndGetWithZeroOptimization) {
    SparseMatrix<double> m(3, 3);

    m.set(0, 0, 5.0);
    m.set(1, 2, -3.5);

    EXPECT_DOUBLE_EQ(m.get(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(m.get(1, 2), -3.5);
    EXPECT_EQ(m.get_non_zero_count(), 2);

    EXPECT_DOUBLE_EQ(m.get(2, 2), 0.0);
    EXPECT_DOUBLE_EQ(m.get(0, 1), 0.0);

    m.set(0, 0, 0.0);
    EXPECT_DOUBLE_EQ(m.get(0, 0), 0.0);
    EXPECT_EQ(m.get_non_zero_count(), 1); 
}

TEST(SparseMatrixTest, OutOfBounds) {
    SparseMatrix<double> m(2, 2);

    EXPECT_THROW(m.get(-1, 0), std::invalid_argument);
    EXPECT_THROW(m.get(0, 2), std::invalid_argument); 
    EXPECT_THROW(m.set(2, 0, 5.0), std::invalid_argument);
}

TEST(SparseMatrixTest, Norm) {
    SparseMatrix<double> m(3, 3);
    m.set(0, 0, 3.0);
    m.set(1, 1, -4.0);
    // norm = sqrt(3^2 + (-4)^2) = sqrt(9 + 16) = 5
    EXPECT_DOUBLE_EQ(m.norm(), 5.0);
}

TEST(SparseMatrixTest, ScalarMultiplication) {
    SparseMatrix<double> m(2, 2);
    m.set(0, 1, 2.5);
    m.set(1, 0, -1.0);

    IMatrix<double>* res = m.mult_scalar(2.0);
    auto* sparse_res = dynamic_cast<SparseMatrix<double>*>(res);

    ASSERT_NE(sparse_res, nullptr);
    EXPECT_DOUBLE_EQ(sparse_res->get(0, 1), 5.0);
    EXPECT_DOUBLE_EQ(sparse_res->get(1, 0), -2.0);
    EXPECT_EQ(sparse_res->get_non_zero_count(), 2);

    IMatrix<double>* res_zero = m.mult_scalar(0.0);
    auto* sparse_zero = dynamic_cast<SparseMatrix<double>*>(res_zero);
    EXPECT_EQ(sparse_zero->get_non_zero_count(), 0);

    delete res;
    delete res_zero;
}

TEST(SparseMatrixTest, Addition) {
    SparseMatrix<double> m1(2, 2);
    m1.set(0, 0, 1.0);
    m1.set(0, 1, 2.0);

    SparseMatrix<double> m2(2, 2);
    m2.set(0, 1, -2.0); 
    m2.set(1, 1, 5.0);

    IMatrix<double>* res = m1.add(m2);
    auto* sparse_res = dynamic_cast<SparseMatrix<double>*>(res);

    ASSERT_NE(sparse_res, nullptr);
    EXPECT_DOUBLE_EQ(sparse_res->get(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(sparse_res->get(0, 1), 0.0); // 2.0 + (-2.0) = 0.0
    EXPECT_DOUBLE_EQ(sparse_res->get(1, 1), 5.0);

    EXPECT_EQ(sparse_res->get_non_zero_count(), 2); 

    SparseMatrix<double> m3(3, 3);
    EXPECT_THROW(m1.add(m3), std::invalid_argument);

    delete res;
}

TEST(SparseMatrixTest, MatrixMultiplication) {
    // m1 size 2x3
    SparseMatrix<double> m1(2, 3);
    m1.set(0, 1, 2.0);
    m1.set(1, 2, 3.0);

    // m2 size 3x2
    SparseMatrix<double> m2(3, 2);
    m2.set(1, 0, 4.0);
    m2.set(2, 1, 5.0);

    // awaiting for reault (2x2):
    // (0,0) = m1(0,1)*m2(1,0) = 2*4 = 8
    // (1,1) = m1(1,2)*m2(2,1) = 3*5 = 15
    IMatrix<double>* res = m1.mult(m2);
    auto* sparse_res = dynamic_cast<SparseMatrix<double>*>(res);

    ASSERT_NE(sparse_res, nullptr);
    EXPECT_EQ(sparse_res->get_rows(), 2);
    EXPECT_EQ(sparse_res->get_cols(), 2);
    EXPECT_DOUBLE_EQ(sparse_res->get(0, 0), 8.0);
    EXPECT_DOUBLE_EQ(sparse_res->get(1, 1), 15.0);
    EXPECT_DOUBLE_EQ(sparse_res->get(0, 1), 0.0);
    EXPECT_EQ(sparse_res->get_non_zero_count(), 2);

    delete res;
}

TEST(SparseMatrixTest, ComplexNumbersSupport) {
    SparseMatrix<Complex<double>> m(2, 2);
    
    m.set(0, 0, Complex<double>(3.0, 4.0)); 
    
    EXPECT_DOUBLE_EQ(m.get(0, 0).re, 3.0);
    EXPECT_DOUBLE_EQ(m.get(0, 0).im, 4.0);
    
    EXPECT_DOUBLE_EQ(m.get(1, 1).re, 0.0);
    EXPECT_DOUBLE_EQ(m.get(1, 1).im, 0.0);

    EXPECT_DOUBLE_EQ(m.norm(), 5.0);
}