#ifndef __LUP_HPP__
#define __LUP_HPP__

#include "./matrix/matrix.hpp"

template <typename T>
class LUP
{
public:
    typedef Matrix<T> matrix_type;
    typedef typename matrix_type::floating_point floating_point;
    typedef typename matrix_type::value_type value_type;
    typedef typename matrix_type::reference reference;
    typedef typename matrix_type::const_reference const_reference;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::difference_type difference_type;

private:
    static constexpr const floating_point value_epsilon = matrix_type::value_epsilon;
    static constexpr const floating_point value_round_error = matrix_type::value_round_error;
    static constexpr const size_type size_min = matrix_type::size_min;
    static constexpr const size_type size_max = matrix_type::size_max;

    matrix_type l;
    matrix_type u;
    matrix_type p;
    bool det;

    size_type find_max(size_type) const;
    void column_to_null(size_type);

    void forward(const matrix_type &, size_type, matrix_type &, size_type) const;
    void back(const matrix_type &, size_type, matrix_type &, size_type) const;

public:
    LUP() = delete;
    LUP(const matrix_type &);
    LUP(const LUP &) = default;
    LUP(LUP &&) = default;
    LUP &operator =(const LUP &) = default;
    LUP &operator =(LUP &&) = default;
    ~LUP() = default;

    const matrix_type &get_u() const;
    const matrix_type &get_l() const;
    const matrix_type &get_p() const;

    value_type determinant() const;
    matrix_type solution(const matrix_type &) const;
    matrix_type invertible() const;
};

template <typename T>
LUP<T>::LUP(const LUP<T>::matrix_type &m) : l(), u(), p(), det()
{
    if (m.size1() != m.size2() || !m.size1())
    {
        throw std::domain_error("Non-square matrix");
    }

    const size_type size = m.size1(), iteration_cnt = m.size1() - 1u;
    u = m;
    l.zero(size, size);
    p.identity(size);

    for (size_type i = 0u; i < iteration_cnt; ++i)
    {
        size_type max = find_max(i);
        if (max != i)
        {
            l.row_switching(max, i);
            u.row_switching(max, i);
            p.row_switching(max, i);
            det ^= true;
        }

        column_to_null(i);
    }

    for (size_type i = 0u; i < size; ++i)
    {
        l(i, i) = 1.0;
    }
}

template <typename T>
typename LUP<T>::size_type LUP<T>::find_max(const LUP<T>::size_type idx) const
{
    const size_type size = u.size1();
    size_type max = size_max;
    for (size_type i = idx; i < size; ++i)
    {
        if (std::isgreaterequal(std::fabs(u(i, idx)), value_epsilon) &&
            (max == size_max ||
                std::isless(std::fabs(u(max, idx)), std::fabs(u(i, idx)))))
        {
            max = i;
        }
    }
    if (max == size_max)
    {
        throw std::underflow_error("Non-invertible matrix");
    }

    return max;
}

template <typename T>
void LUP<T>::column_to_null(const LUP<T>::size_type idx)
{
    const size_type size = u.size1();
    for (size_type i = idx + 1u; i < size; ++i)
    {
        const value_type mu = -u(i, idx) / u(idx, idx);
        if (std::fabs(mu) >= value_epsilon)
        {
            l(i, idx) = -mu;
            u.row_addition(i, mu, idx);
        }
    }
}

template <typename T>
void LUP<T>::forward(const typename LUP<T>::matrix_type &b,
    const typename LUP<T>::size_type b_idx,
    typename LUP<T>::matrix_type &x,
    const typename LUP<T>::size_type x_idx) const
{
    const difference_type size = l.size1();
    for (difference_type i = 0; i < size; ++i)
    {
        x(i, x_idx) = b(i, b_idx);
        for (difference_type j = 0; j < i; ++j)
        {
            x(i, x_idx) -= l(i, j) * x(j, x_idx);
        }
    }
}

template <typename T>
void LUP<T>::back(const typename LUP<T>::matrix_type &b,
    const typename LUP<T>::size_type b_idx,
    typename LUP<T>::matrix_type &x,
    const typename LUP<T>::size_type x_idx) const
{
    const difference_type size = u.size1();
    for (difference_type i = size - 1; i >= 0; --i)
    {
        x(i, x_idx) = b(i, b_idx);
        for (difference_type j = size - 1; j > i; --j)
        {
            x(i, x_idx) -= u(i, j) * x(j, x_idx);
        }
        x(i, x_idx) /= u(i, i);
    }
}

template <typename T>
const typename LUP<T>::matrix_type &LUP<T>::get_u() const
{
    return u;
}

template <typename T>
const typename LUP<T>::matrix_type &LUP<T>::get_l() const
{
    return l;
}

template <typename T>
const typename LUP<T>::matrix_type &LUP<T>::get_p() const
{
    return p;
}

template <typename T>
typename LUP<T>::value_type LUP<T>::determinant() const
{
    const size_type size = u.size1();
    value_type result = det ? -1.0 : 1.0;
    for (size_type i = 0u; i < size; ++i)
    {
        result *= u(i, i);
    }

    return result;
}

template <typename T>
typename LUP<T>::matrix_type LUP<T>::solution(const typename LUP<T>::matrix_type &b) const
{
    if (l.size1() != b.size1() || b.size2() != 1u)
    {
        throw std::domain_error("Can't be solved");
    }

    const Matrix<T> p_mul_b = p * b;
    Matrix<T> x(u.size1(), 1u);
    forward(p_mul_b, 0u, x, 0u);
    back(x, 0u, x, 0u);

    return x;
}

template <typename T>
typename LUP<T>::matrix_type LUP<T>::invertible() const
{
    Matrix<T> x(l.size1(), u.size2());
    const size_type size = l.size1();
    for (size_type i = 0u; i < size; ++i)
    {
        forward(p, i, x, i);
    }
    for (size_type i = 0u; i < size; ++i)
    {
        back(x, i, x, i);
    }

    return x;
}

#endif
