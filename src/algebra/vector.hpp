#pragma once

#include "../lab2/memory/dynamic_array.hpp"
#include "../lab2/types/exceptions.hpp"
#include "../types/math.hpp"

template <typename T>
class Vector {
private:
    DynamicArray<T> m_data;

public:
    // constructors
    Vector();
    explicit Vector(int size);
    Vector(int size, const T& fill_value);
    Vector(T* items, int count);

    // The Rule of Zero (Works perfectly, as there are no primitive size fields)
    ~Vector() = default;
    Vector(const Vector<T>& other) = default;
    Vector(Vector<T>&& other) noexcept = default;
    Vector<T>& operator=(const Vector<T>& other) = default;
    Vector<T>& operator=(Vector<T>&& other) noexcept = default;

    // getters
    [[nodiscard]] int get_size() const;
    [[nodiscard]] const T& get(int index) const;
    void set(int index, const T& value);

    // Overloading the [] operator
    const T& operator[](int index) const; // reading only: val = v[0]

    //TODO: += *= штуки добавить для оптимизации
    // math operations
    Vector<T> operator+(const Vector<T>& other) const;
    Vector<T> operator-(const Vector<T>& other) const;
    Vector<T> operator*(const T& scalar) const;

    // Scalar Product of Vectors (Dot Product)
    T dot(const Vector<T>& other) const;

    // Euclidean norm (vector length)
    double norm() const;
};
//TODO: можно использовать % для скалярного произведения 
//TODO: добавить векторное произведение? 
// Global operator for a scalar on the left (scalar * vector)
template <typename T>
Vector<T> operator*(const T& scalar, const Vector<T>& v);

#include "vector.tpp"