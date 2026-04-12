#pragma once
#include "complex.hpp"

namespace custom_math {
    // 1. Base template for standard numeric types (int, double, float)
    // T{} is a safe zero for any type
    template <typename T>
    constexpr T abs(const T& value) {
        return (value < T{}) ? -value : value;
    }

    // 2. Overload specifically for complex numbers
    // The compiler will automatically select this function if a Complex<T> is passed.
    template <typename T>
    constexpr double abs(const Complex<T>& c) {
        return c.abs();
    }

}