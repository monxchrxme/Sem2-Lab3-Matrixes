#pragma once

/**
 * @brief abstract interface for matrix
 * * requirements to T type:
 * - T{} — have to return the additive neutral element (zero)
 * - T{1} — have to return the multiply neutral element (one)
 * - The +, -, *, and modulo (or squaring) operators must be overload
 * for calculating the norm
 */
template <typename T>
class IMatrix {
public:
    virtual ~IMatrix() {}

    // getters
    virtual int get_rows() const = 0;
    virtual int get_cols() const = 0;

    // method get for polymorphic reading
    virtual const T& get(int i, int j) const = 0;

    // access operator (only reading)
    const T& operator()(int i, int j) const {
        return get(i, j);
    }

    // operations
    virtual void set(int i, int j, const T &value) = 0;

    // math operations
    virtual IMatrix<T>* add(const IMatrix<T> &other) const = 0;    // returns a heap-allocated object, must be deleted by user
    virtual IMatrix<T>* mult_scalar(const T &value) const = 0;     // returns a heap-allocated object, must be deleted by user
    virtual IMatrix<T>* mult(const IMatrix<T> &other) const = 0;   // returns a heap-allocated object, must be deleted by user

    // Frobenius's norm calculating
    virtual double norm() const = 0;

    // method for deep copy
    virtual IMatrix<T>* clone() const = 0;
};