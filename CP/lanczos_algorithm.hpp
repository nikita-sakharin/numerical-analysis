#ifndef __LANCZOS_ALGORITHM_HPP__
#define __LANCZOS_ALGORITHM_HPP__

#include <cstddef>

#include <exception>
#include <stdexcept>
#include <type_traits>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric;

template<typename T>
using SparseMatrix = ublas::compressed_matrix<T>;

template<typename T>
using SparseVector = ublas::compressed_vector<T>;

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
SparseVector<T> lanczos_algorithm(const SparseMatrix<T> &a, const SparseVector<T> &b);

template <typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
static void qr_decomposition(const SparseMatrix<T> &,
    SparseMatrix<T> &, SparseMatrix<T> &);

template <typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
static void householder_matrix(const SparseVector<T> &, SparseMatrix<T> &);

template <typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
static void get_v(const SparseMatrix<T> &, std::size_t,
    SparseVector<T> &v);

template <typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
static void tridiagonal_eigenvalue(const SparseMatrix<T> &, const T &, const SparseMatrix<T> &, SparseVector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
SparseVector<T> lanczos_algorithm(const SparseMatrix<T> &a, const SparseVector<T> &b)
{
    static constexpr T
        EPS = std::pow(static_cast<T>(10.), -std::numeric_limits<T>::digits10);
    const std::size_t size = a.size1();
    if (size != a.size2() || size != b.size())
    {
        throw std::logic_error("Matrix A must be square");
    }

    SparseMatrix<T> identity(size, size);
    for (std::size_t i = 0; i < size; ++i)
    {
        identity(i, i) = 1.;
    }

    SparseMatrix<T> q(size, size);
    ublas::column(q, 0) = b / ublas::norm_2(b);
    SparseVector<T> z;
    T alpha, beta = 0.;
    for (std::size_t i = 0; i < size; ++i)
    {
        z = ublas::prod(a, ublas::column(q, i));
        alpha = ublas::inner_prod(ublas::column(q, i), z);

        z -= alpha * ublas::column(q, i);
        if (i)
        {
            z -= beta * ublas::column(q, i - 1);
        }

        beta = ublas::norm_2(z);
        if (beta < EPS)
        {
            break;
        }
        if (i + 1 < size)
        {
            ublas::column(q, i + 1) = z / beta;
        }
    }

    SparseMatrix<T> tridiag = ublas::prod(ublas::trans(q), a);
    tridiag = ublas::prod(tridiag, q);
    SparseVector<T> result;
    tridiagonal_eigenvalue(tridiag, EPS, identity, result);

    return result;
}

template <typename T, typename>
static void tridiagonal_eigenvalue(const SparseMatrix<T> &tridiag, const T &eps,
    const SparseMatrix<T> &identity, SparseVector<T> &result)
{
    static constexpr std::size_t LIMIT = 2;
    SparseMatrix<T> t, q, r;
    const std::size_t size = tridiag.size1();
    result.resize(size);
    for (std::size_t m = size; m > LIMIT; --m)
    {
        const std::size_t i = m - 1;
        t = ublas::subrange(tridiag, 0, m, 0, m);
        while (std::fabs(t(i, i - 2)) >= eps * std::fabs(t(i, i)))
        {
            const T mu = t(i, i);
            t -= mu * ublas::subrange(identity, 0, m, 0, m);
            qr_decomposition(t, q, r);
            axpy_prod(r, q, t);
            t += mu * ublas::subrange(identity, 0, m, 0, m);
        }
        const T value = t(i, i);
        if (std::fabs(value) >= eps)
        {
            result(i) = value;
        }
    }
    for (std::size_t i = 0; i < LIMIT; ++i)
    {
        const T value = t(i, i);
        if (std::fabs(value) >= eps)
        {
            result(i) = value;
        }
    }
}

template <typename T, typename>
static void qr_decomposition(const SparseMatrix<T> &a, SparseMatrix<T> &q, SparseMatrix<T> &r)
{
    const std::size_t size = a.size1();
    if (!size || size != a.size2())
    {
        throw std::logic_error("Non-square matrix");
    }

    r = a;
    q = ublas::identity_matrix<T>(size);

    SparseVector<T> v(size);
    SparseMatrix<T> h_k;
    for (std::size_t i = 0u; i < size - 1u; ++i)
    {
        get_v(r, i, v);
        householder_matrix(v, h_k);
        q = ublas::prod(q, h_k);
        r = ublas::prod(h_k, r);
    }
}

template <typename T, typename>
static void householder_matrix(const SparseVector<T> &v,
    SparseMatrix<T> &h_k)
{
    h_k = ublas::outer_prod(v, v);
    const T v_multiplies_v_t = -2. / ublas::inner_prod(v, v);
    h_k *= v_multiplies_v_t;
    h_k += ublas::identity_matrix<T>(v.size());
}

template <typename T, typename>
static void get_v(const SparseMatrix<T> &r, const std::size_t idx,
    SparseVector<T> &v)
{
    const std::size_t size = r.size1();
    if (idx)
    {
        v.erase_element(idx - 1);
    }
    ublas::subrange(v, idx, size) = ublas::subrange(ublas::column(r, idx), idx, size);
    v(idx) += std::copysign(ublas::norm_2(v), r(idx, idx));
}

#endif
