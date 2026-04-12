#pragma once

#include <cmath>
#include <stdexcept>

template <typename T>
struct Complex {
    T re;
    T im;

    Complex(T re = T{}, T im = T{}) : re{re}, im{im} {}

/***********************************************************************
 * OPERATIONS
 ***********************************************************************/

    // unary minus
    Complex operator-() const {
        return Complex<T>(-re, -im);
    }

    // binary operations
    Complex operator+(const Complex& o) const {
        return Complex<T>(re + o.re, im + o.im);
    }

    Complex operator-(const Complex& o) const {
        return Complex<T>(re - o.re, im - o.im);
    }

    Complex operator*(const Complex& o) const {
        return Complex<T>(re * o.re - im * o.im, re * o.im + im * o.re);
    }

    Complex operator/(const Complex& o) const {
        T div = o.re * o.re + o.im * o.im;
        if (div == T{}) {
            throw std::domain_error("Complex::operator/: division by zero");
        }
        return Complex<T>((re * o.re + im * o.im) / div, (im * o.re - re * o.im) / div);
    }

    Complex& operator+=(const Complex& o) {
        re += o.re;
        im += o.im;
        return *this;
    }

    Complex& operator-=(const Complex& o) {
        re -= o.re;
        im -= o.im;
        return *this;
    }

    Complex& operator*=(const Complex& o) {
        T new_re = re * o.re - im * o.im;
        T new_im = re * o.im + im * o.re;
        re = new_re;
        im = new_im;
        return *this;
    }

    Complex& operator/=(const Complex& o) {
        *this = *this / o;
        return *this;
    }

    // comparison operators
    bool operator==(const Complex& o) const {
        return re == o.re && im == o.im;
    }

    bool operator!=(const Complex& o) const {
        return !(*this == o);
    }

    // math funcs
    T norm_sqr() const {
        return re * re + im * im;
    }

    // returns a double for precise root calculation
    double abs() const {
        return std::sqrt(static_cast<double>(norm_sqr()));
    }
};

// global operators (scalar on the left)
template <typename T>
Complex<T> operator+(const T& scalar, const Complex<T>& c) {
    return Complex<T>(scalar) + c;
}

template <typename T>
Complex<T> operator-(const T& scalar, const Complex<T>& c) {
    return Complex<T>(scalar) - c;
}

template <typename T>
Complex<T> operator*(const T& scalar, const Complex<T>& c) {
    return Complex<T>(scalar) * c;
}

template <typename T>
Complex<T> operator/(const T& scalar, const Complex<T>& c) {
    return Complex<T>(scalar) / c;
}

// testing utility (Floating-Point Comparison)
template <typename T>
bool approx_equal(const Complex<T>& a, const Complex<T>& b, double eps = 1e-9) {
    return std::abs(static_cast<double>(a.re - b.re)) < eps &&
           std::abs(static_cast<double>(a.im - b.im)) < eps;
}