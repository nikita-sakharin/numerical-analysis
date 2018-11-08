#ifndef __LUP_HPP__
#define __LUP_HPP__

#include <cmath>
#include <cstddef>

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <type_traits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
void row_switching(ublas::matrix<T> &,
    const std::size_t, const std::size_t) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
void row_multiplication(ublas::matrix<T> &,
    const T, const std::size_t) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
void row_addition(ublas::matrix<T> &,
    const std::size_t, const T, const std::size_t) noexcept;

template <typename T>
class LUP
{
public:
    typedef typename std::enable_if<std::is_floating_point<T>::value, T>::type
        floating_point;

    typedef ublas::vector<T> vector_type;
    typedef ublas::matrix<T> matrix_type;

    typedef typename matrix_type::value_type value_type;
    typedef typename matrix_type::reference reference;
    typedef typename matrix_type::const_reference const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

private:
    static constexpr floating_point
        value_epsilon = std::numeric_limits<floating_point>::epsilon(),
        value_round_error = std::numeric_limits<floating_point>::round_error();
    static constexpr std::size_t
        size_min = std::numeric_limits<std::size_t>::min(),
        size_max = std::numeric_limits<std::size_t>::max();
    static constexpr std::ptrdiff_t
        difference_min = std::numeric_limits<std::ptrdiff_t>::min(),
        difference_max = std::numeric_limits<std::ptrdiff_t>::max();

    matrix_type l;
    matrix_type u;
    matrix_type p;
    bool det;

    std::size_t find_max(std::size_t) const;
    void column_to_null(std::size_t) noexcept;

    void forward(const matrix_type &, std::size_t,
        matrix_type &, std::size_t) const noexcept;
    void back(const matrix_type &, std::size_t,
        matrix_type &, std::size_t) const noexcept;

public:
    LUP() = delete;
    LUP(const matrix_type &);
    LUP(const LUP &) = default;
    LUP(LUP &&) = default;
    LUP &operator =(const LUP &) = default;
    LUP &operator =(LUP &&) = default;
    ~LUP() = default;

    const matrix_type &get_u() const noexcept;
    const matrix_type &get_l() const noexcept;
    const matrix_type &get_p() const noexcept;

    value_type determinant() const noexcept;
    vector_type solution(const vector_type &) const;
    matrix_type invertible() const noexcept;
};

template <typename T>
LUP<T>::LUP(const LUP<T>::matrix_type &m) : l(), u(), p(), det()
{
    if (m.size1() != m.size2() || !m.size1())
    {
        throw std::domain_error("Non-square matrix");
    }

    const std::size_t size = m.size1(), iteration_cnt = m.size1() - 1;
    u = m;
    l = ublas::zero_matrix<T>(size, size);
    p = ublas::identity_matrix<T>(size);

    for (std::size_t i = 0; i < iteration_cnt; ++i)
    {
        std::size_t max = find_max(i);
        if (max != i)
        {
            row_switching(l, max, i);
            row_switching(u, max, i);
            row_switching(p, max, i);
            det ^= true;
        }

        column_to_null(i);
    }

    for (std::size_t i = 0; i < size; ++i)
    {
        l(i, i) = 1.0;
    }
}

template <typename T>
std::size_t LUP<T>::find_max(const std::size_t idx) const
{
    const std::size_t size = u.size1();
    std::size_t max = size_max;
    for (std::size_t i = idx; i < size; ++i)
    {
        if (std::isgreaterequal(std::abs(u(i, idx)), value_epsilon) &&
            (max == size_max ||
                std::isless(std::abs(u(max, idx)), std::abs(u(i, idx)))))
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
void LUP<T>::column_to_null(const std::size_t idx) noexcept
{
    const std::size_t size = u.size1();
    for (std::size_t i = idx + 1; i < size; ++i)
    {
        const value_type mu = -u(i, idx) / u(idx, idx);
        if (std::abs(mu) >= value_epsilon)
        {
            l(i, idx) = -mu;
            row_addition(u, i, mu, idx);
        }
    }
}

template <typename T>
void LUP<T>::forward(const LUP<T>::matrix_type &b,
    const std::size_t b_idx,
    LUP<T>::matrix_type &x,
    const std::size_t x_idx) const noexcept
{
    const std::ptrdiff_t size = l.size1();
    for (std::ptrdiff_t i = 0; i < size; ++i)
    {
        x(i, x_idx) = b(i, b_idx);
        for (std::ptrdiff_t j = 0; j < i; ++j)
        {
            x(i, x_idx) -= l(i, j) * x(j, x_idx);
        }
    }
}

template <typename T>
void LUP<T>::back(const LUP<T>::matrix_type &b,
    const std::size_t b_idx,
    LUP<T>::matrix_type &x,
    const std::size_t x_idx) const noexcept
{
    const std::ptrdiff_t size = u.size1();
    for (std::ptrdiff_t i = size - 1; i >= 0; --i)
    {
        x(i, x_idx) = b(i, b_idx);
        for (std::ptrdiff_t j = size - 1; j > i; --j)
        {
            x(i, x_idx) -= u(i, j) * x(j, x_idx);
        }
        x(i, x_idx) /= u(i, i);
    }
}

template <typename T>
const typename LUP<T>::matrix_type &LUP<T>::get_u() const noexcept
{
    return u;
}

template <typename T>
const typename LUP<T>::matrix_type &LUP<T>::get_l() const noexcept
{
    return l;
}

template <typename T>
const typename LUP<T>::matrix_type &LUP<T>::get_p() const noexcept
{
    return p;
}

template <typename T>
typename LUP<T>::value_type LUP<T>::determinant() const noexcept
{
    const std::size_t size = u.size1();
    value_type result = det ? -1.0 : 1.0;
    for (std::size_t i = 0; i < size; ++i)
    {
        result *= u(i, i);
    }

    return result;
}

template <typename T>
typename LUP<T>::vector_type LUP<T>::solution(const LUP<T>::vector_type &b) const
{
    if (l.size1() != b.size())
    {
        throw std::domain_error("Can't be solved");
    }

    matrix_type p_mul_b(u.size1(), 1), x(u.size1(), 1);
    column(p_mul_b, 0) = prod(p, b);
    forward(p_mul_b, 0, x, 0);
    back(x, 0, x, 0);

    return column(x, 0);
}

template <typename T>
typename LUP<T>::matrix_type LUP<T>::invertible() const noexcept
{
    matrix_type x(l.size1(), u.size2());
    const std::size_t size = l.size1();
    for (std::size_t i = 0; i < size; ++i)
    {
        forward(p, i, x, i);
    }
    for (std::size_t i = 0; i < size; ++i)
    {
        back(x, i, x, i);
    }

    return x;
}

template <typename T, typename>
void row_switching(ublas::matrix<T> &m,
    const std::size_t i, const std::size_t j) noexcept
{
    const std::size_t size = m.size2();
    for (std::size_t k = 0; k < size; ++k)
    {
        std::swap(m(i, k), m(j, k));
    }
}

template <typename T, typename>
void row_multiplication(ublas::matrix<T> &m,
    const T alpha, const std::size_t i) noexcept
{
    const std::size_t size = m.size2();
    for (std::size_t k = 0; k < size; ++k)
    {
        m(i, k) *= alpha;
    }
}

template <typename T, typename>
void row_addition(ublas::matrix<T> &m,
    const std::size_t i, const T alpha, const std::size_t j) noexcept
{
    const std::size_t size = m.size2();
    for (std::size_t k = 0; k < size; ++k)
    {
        m(i, k) = std::fma(alpha, m(j, k), m(i, k));
    }
}

#endif
