#include <gtest/gtest.h>
#include "../src/types/complex.hpp"
#include <stdexcept>

TEST(ComplexTest, Constructors) {
Complex<double> c_default;
EXPECT_DOUBLE_EQ(c_default.re, 0.0);
EXPECT_DOUBLE_EQ(c_default.im, 0.0);

Complex<double> c_scalar(5.5);
EXPECT_DOUBLE_EQ(c_scalar.re, 5.5);
EXPECT_DOUBLE_EQ(c_scalar.im, 0.0);

Complex<double> c_full(1.2, 3.4);
EXPECT_DOUBLE_EQ(c_full.re, 1.2);
EXPECT_DOUBLE_EQ(c_full.im, 3.4);
}

TEST(ComplexTest, EqualityAndApproximation) {
Complex<double> c1(1.000000001, 2.0);
Complex<double> c2(1.0, 2.0);

EXPECT_NE(c1, c2);
EXPECT_TRUE(approx_equal(c1, c2, 1e-6));
}

TEST(ComplexTest, AdditionAndSubtraction) {
Complex<double> c1(2.0, 3.0);
Complex<double> c2(4.0, -1.0);

Complex<double> sum = c1 + c2;
EXPECT_DOUBLE_EQ(sum.re, 6.0);
EXPECT_DOUBLE_EQ(sum.im, 2.0);

Complex<double> diff = c1 - c2;
EXPECT_DOUBLE_EQ(diff.re, -2.0);
EXPECT_DOUBLE_EQ(diff.im, 4.0);

Complex<double> neg = -c1;
EXPECT_DOUBLE_EQ(neg.re, -2.0);
EXPECT_DOUBLE_EQ(neg.im, -3.0);
}

TEST(ComplexTest, ScalarArithmetic) {
Complex<double> c(2.0, 3.0);

Complex<double> add_right = c + 5.0;
Complex<double> add_left = 5.0 + c;
EXPECT_EQ(add_right, Complex<double>(7.0, 3.0));
EXPECT_EQ(add_left, add_right);

Complex<double> mult_right = c * 2.0;
Complex<double> mult_left = 2.0 * c;
EXPECT_EQ(mult_right, Complex<double>(4.0, 6.0));
EXPECT_EQ(mult_left, mult_right);

Complex<double> sub_left = 10.0 - c; // (10+0i) - (2+3i) = 8 - 3i
EXPECT_EQ(sub_left, Complex<double>(8.0, -3.0));
}

TEST(ComplexTest, Multiplication) {
Complex<double> c1(1.0, 2.0);
Complex<double> c2(3.0, 4.0);
// (1+2i)(3+4i) = 3 + 4i + 6i - 8 = -5 + 10i
Complex<double> mult = c1 * c2;
EXPECT_DOUBLE_EQ(mult.re, -5.0);
EXPECT_DOUBLE_EQ(mult.im, 10.0);
}

TEST(ComplexTest, DivisionAndEdgeCases) {
Complex<double> c1(1.0, 2.0);
Complex<double> c2(3.0, 4.0);
// (1+2i)/(3+4i) = (3 + 8) / 25 + i(6 - 4) / 25 = 11/25 + 2/25i
Complex<double> div = c1 / c2;
EXPECT_DOUBLE_EQ(div.re, 0.44);
EXPECT_DOUBLE_EQ(div.im, 0.08);
Complex<double> zero(0.0, 0.0);
EXPECT_THROW(c1 / zero, std::domain_error);
EXPECT_THROW(c1 / 0.0, std::domain_error);
}

TEST(ComplexTest, CompoundAssignments) {
Complex<double> c(2.0, 3.0);
c += Complex<double>(1.0, 1.0);
EXPECT_EQ(c, Complex<double>(3.0, 4.0));

c *= 2.0;
EXPECT_EQ(c, Complex<double>(6.0, 8.0));
}

TEST(ComplexTest, MathFunctions) {
Complex<double> c(3.0, 4.0);
EXPECT_DOUBLE_EQ(c.norm_sqr(), 25.0);
EXPECT_DOUBLE_EQ(c.abs(), 5.0);
}

TEST(ComplexTest, IntegerTemplateType) {
Complex<int> c1(3, 4);
Complex<int> c2(1, 2);

Complex<int> sum = c1 + c2;
EXPECT_EQ(sum.re, 4);
EXPECT_EQ(sum.im, 6);

EXPECT_DOUBLE_EQ(c1.abs(), 5.0);
}