#ifndef __THOMAS_ALGORITHM_HPP__
#define __THOMAS_ALGORITHM_HPP__

#include "./matrix/matrix.hpp"

template <typename T>
Matrix<T> thomas_algorithm(const Matrix<T> &a, const Matrix<T> &b, const Matrix<T> &c,
    const Matrix<T> &d)
{
    if (a.size1() != b.size1() || b.size1() != c.size1() || c.size1() != d.size1() ||
        a.size1() != 1u ||
        a.size2() != b.size2() || b.size2() != c.size2() || c.size2() != d.size2() ||
        !d.size2())
    {
        throw std::domain_error("Can't be solved");
    }

    const typename Matrix<T>::size_type size = d.size2();
    Matrix<T> p(1u, size), q(1u, size), x(size, 1u);

    p(0u) = -c(0u) / b(0u);
    q(0u) = d(0u) / b(0u);
    for (typename Matrix<T>::size_type i = 1u; i < size; ++i)
    {
        p(i) = -c(i) / (b(i) + a(i) * p(i - 1u));
        q(i) = (d(i) - a(i) * q(i - 1u)) / (b(i) + a(i) * p(i - 1u));
    }

    x(size - 1u) = q(size - 1u);
    for (typename Matrix<T>::difference_type i = size - 2u; i >= 0; --i)
    {
        x(i) = p(i) * x(i + 1u) + q(i);
    }

    return x;
}

#endif
