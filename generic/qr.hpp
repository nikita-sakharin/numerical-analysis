#ifndef __QR_HPP__
#define __QR_HPP__

#include <cmath>
#include <cstddef>

#include <complex>
#include <exception>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<std::complex<T>> qr(const ublas::matrix<T> &, const T &);

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static inline void arg_vilid(std::size_t, std::size_t, const T &);

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void find_decomposition(const ublas::matrix<T> &,
    ublas::matrix<T> &, ublas::matrix<T> &) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void householder_matrix(const ublas::vector<T> &, ublas::matrix<T> &) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void get_v(const ublas::matrix<T> &, std::size_t, ublas::vector<T> &) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static T criteria(const ublas::matrix<T> &) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void find_root(const ublas::matrix<T> &,
    ublas::vector<std::complex<T>> &, const T &) noexcept;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static std::pair<std::complex<T>, std::complex<T>>
    quadratic_equation(const T &, const T&, const T &) noexcept;

template <typename T, typename>
ublas::vector<std::complex<T>> qr(const ublas::matrix<T> &a, const T &eps)
{
    arg_vilid<T>(a.size1(), a.size2(), eps);
    ublas::matrix<T> q, r, a_k = a;
    do
    {
        find_decomposition(a_k, q, r);
        ublas::axpy_prod(r, q, a_k);
    } while (criteria(a_k) > eps);

    ublas::vector<std::complex<T>> result;
    find_root(a_k, result, eps);

    return result;
}

template <typename T, typename>
static inline void arg_vilid(const std::size_t m, const std::size_t n,
    const T &eps)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    if (!m || m != n || eps < EPSILON)
    {
        throw std::logic_error("Non-square matrix or eps is too small");
    }
}

template <typename T, typename>
static void find_decomposition(const ublas::matrix<T> &a, ublas::matrix<T> &q, ublas::matrix<T> &r) noexcept
{
    const std::size_t m = a.size1();
    r = a;
    q = ublas::identity_matrix<T>(m);

    ublas::vector<T> v;
    ublas::matrix<T> h_k;
    for (std::size_t i = 0; i < m - 1; ++i)
    {
        get_v(r, i, v);
        householder_matrix(v, h_k);
        q = ublas::prod(q, h_k);
        r = ublas::prod(h_k, r);
    }
}

template <typename T, typename>
static void householder_matrix(const ublas::vector<T> &v,
    ublas::matrix<T> &h_k) noexcept
{
    h_k = ublas::outer_prod(v, v);
    const T v_multiplies_v_t = -2.0 / ublas::inner_prod(v, v);
    h_k *= v_multiplies_v_t;
    h_k += ublas::identity_matrix<T>(v.size());
}

template <typename T, typename>
static void get_v(const ublas::matrix<T> &r, const std::size_t idx, ublas::vector<T> &v) noexcept
{
    const std::size_t m = r.size1();
    v.resize(m);

    T sum = 0.0;
    for (std::size_t i = 0; i < m; ++i)
    {
        if (i < idx)
        {
            v(i) = 0.0;
        }
        else
        {
            v(i) = r(i, idx);
            sum += r(i, idx) * r(i, idx);
        }
    }
    v(idx) += std::copysign(std::sqrt(sum), r(idx, idx));
}

template <typename T, typename>
static T criteria(const ublas::matrix<T> &a_k) noexcept
{
    const std::size_t m = a_k.size1();
    T sum = 0.0;
    for (std::size_t i = 2; i < m; ++i)
    {
        for (std::size_t j = 0; j < m - 2; ++j)
        {
            sum += a_k(i, j) * a_k(i, j);
        }
    }

    return std::sqrt(sum);
}


template <typename T, typename>
static void find_root(const ublas::matrix<T> &a_k,
    ublas::vector<std::complex<T>> &result, const T &eps) noexcept
{
    const std::size_t m = a_k.size1();
    result.resize(m);
    std::pair<std::complex<T>, std::complex<T>> root;
    for (std::size_t i = 0; i < m; ++i)
    {
        if (i < m - 1 && a_k(i + 1, i) >= eps)
        {
            const T a = 1.0, b = -(a_k(i, i) + a_k(i + 1, i + 1)),
                c = a_k(i, i) * a_k(i + 1, i + 1) - a_k(i, i + 1) * a_k(i + 1, i);
            root = quadratic_equation<T>(a, b, c);
            result(i) = root.first;
            ++i;
            result(i) = root.second;
        }
        else
        {
            result(i) = a_k(i, i);
        }
    }
}

template <typename T, typename>
std::pair<std::complex<T>, std::complex<T>>
    quadratic_equation(const T &a, const T &b, const T &c) noexcept
{
    const std::complex<T> d = b * b - T(4.0) * a * c,
        a_multiplies_two = T(2.0) * a,
        x1 = (-b + std::sqrt(d)) / a_multiplies_two,
        x2 = (-b - std::sqrt(d)) / a_multiplies_two;

    return std::make_pair(x1, x2);
}

#endif
