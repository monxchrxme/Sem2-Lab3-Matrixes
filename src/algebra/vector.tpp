#pragma once
#include "vector.hpp"
#include <cmath>
#include <stdexcept>

// constructors
template <typename T>
Vector<T>::Vector() : m_data() {}

template <typename T>
Vector<T>::Vector(int size) : m_data(size) {
    if (size < 0) {
        throw std::invalid_argument("Vector size cannot be negative");
    }
}

template <typename T>
Vector<T>::Vector(int size, const T& fill_value) : m_data(size) {
    if (size < 0) {
        throw std::invalid_argument("Vector size cannot be negative");
    }
    for (int i = 0; i < size; ++i) {
        m_data.set(i, fill_value);
    }
}

template <typename T>
Vector<T>::Vector(T* items, int count) : m_data(items, count) {
    if (count < 0) {
        throw std::invalid_argument("Vector size cannot be negative");
    }
}

// getters

template <typename T>
int Vector<T>::get_size() const {
    return m_data.get_size();
}

template <typename T>
const T& Vector<T>::get(int index) const {
    if (index < 0 || index >= get_size()) {
        throw IndexOutOfRange("Vector::get: Index out of bounds");
    }
    return m_data.get(index);
}

template <typename T>
void Vector<T>::set(int index, const T& value) {
    if (index < 0 || index >= get_size()) {
        throw IndexOutOfRange("Vector::set: Index out of bounds");
    }
    m_data.set(index, value);
}

template <typename T>
const T& Vector<T>::operator[](int index) const {
    return this->get(index);
}

// math operations

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& other) const {
    if (get_size() != other.get_size()) {
        throw std::invalid_argument("Vector::operator+: Dimensions mismatch");
    }
    Vector<T> result(get_size());
    for (int i = 0; i < get_size(); ++i) {
        result.set(i, this->get(i) + other.get(i));
    }
    return result;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& other) const {
    if (get_size() != other.get_size()) {
        throw std::invalid_argument("Vector::operator-: Dimensions mismatch");
    }
    Vector<T> result(get_size());
    for (int i = 0; i < get_size(); ++i) {
        result.set(i, this->get(i) - other.get(i));
    }
    return result;
}

template <typename T>
Vector<T> Vector<T>::operator*(const T& scalar) const {
    Vector<T> result(get_size());
    for (int i = 0; i < get_size(); ++i) {
        result.set(i, this->get(i) * scalar);
    }
    return result;
}

template <typename T>
T Vector<T>::dot(const Vector<T>& other) const {
    if (get_size() != other.get_size()) {
        throw std::invalid_argument("Vector::dot: Dimensions mismatch");
    }
    T sum = T{};
    for (int i = 0; i < get_size(); ++i) {
        sum = sum + (this->get(i) * other.get(i));
    }
    return sum;
}

template <typename T>
Vector<T> Vector<T>::cross(const Vector<T>& other) const {
    if (this->get_size() != 3 || other.get_size() != 3) {
        throw std::invalid_argument("Vector::cross: Cross product is strictly defined only for 3D vectors.");
    }

    Vector<T> result(3);
    
    // Calculating coords (indexies 0=X, 1=Y, 2=Z)
    T cx = this->get(1) * other.get(2) - this->get(2) * other.get(1);
    T cy = this->get(2) * other.get(0) - this->get(0) * other.get(2);
    T cz = this->get(0) * other.get(1) - this->get(1) * other.get(0);

    result.set(0, cx);
    result.set(1, cy);
    result.set(2, cz);

    return result;
}

template <typename T>
double Vector<T>::norm() const {
    double sum = 0;
    for (int i = 0; i < get_size(); ++i) {
        double a = static_cast<double>(custom_math::abs(m_data.get(i)));
        sum += a * a;
    }
    return std::sqrt(sum);
}

// Global Operations

template <typename T>
Vector<T> operator*(const T& scalar, const Vector<T>& v) {
    return v * scalar;
}