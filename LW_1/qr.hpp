#ifndef __QR_HPP__
#define __QR_HPP__

#include "./matrix/matrix_complex.hpp"

template <typename T>
static inline void arg_vilid(typename Matrix<T>::size_type, typename Matrix<T>::size_type,
    typename Matrix<T>::floating_point);

template <typename T>
static void find_decomposition(const Matrix<T> &, Matrix<T> &, Matrix<T> &);

template <typename T>
static Matrix<T> householder_matrix(const Matrix<T> &);

template <typename T>
static Matrix<T> get_v(const Matrix<T> &,
    typename Matrix<T>::size_type);

template <typename T>
std::pair<
    std::complex<typename Matrix<T>::floating_point>,
    std::complex<typename Matrix<T>::floating_point>>
quadratic_equation(
    typename Matrix<T>::floating_point,
    typename Matrix<T>::floating_point,
    typename Matrix<T>::floating_point);

template <typename T>
static bool criteria(const Matrix<T> &, typename Matrix<T>::floating_point);

template <typename T>
Matrix<std::complex<typename Matrix<T>::floating_point>> qr(const Matrix<T> &a, typename Matrix<T>::floating_point epsilon)
{
    arg_vilid<T>(a.size1(), a.size2(), epsilon);
    Matrix<T> q, r;
    Matrix<T> a_k = a;
    do
    {
        find_decomposition(a_k, q, r);
        a_k = r * q;
    } while (criteria(a_k, epsilon));
    const typename Matrix<T>::size_type m = a_k.size1();
    Matrix<std::complex<typename Matrix<T>::floating_point>> result(m, 1u);
    std::pair<std::complex<typename Matrix<T>::floating_point>,
        std::complex<typename Matrix<T>::floating_point>> root;
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        if (i + 1u < m && a_k(i + 1u, i) >= epsilon)
        {
            root = quadratic_equation<T>(
                1.0,
                -(a_k(i, i) + a_k(i + 1u, i + 1u)),
                a_k(i, i) * a_k(i + 1u, i + 1u) - a_k(i, i + 1u) * a_k(i + 1u, i));
            result(i) = root.first;
            ++i;
            result(i) = root.second;
        }
        else
        {
            result(i) = a_k(i, i);
        }
    }

    return result;
}

template <typename T>
static inline void arg_vilid(const typename Matrix<T>::size_type a1,
    const typename Matrix<T>::size_type a2,
    const typename Matrix<T>::floating_point epsilon)
{
    if (!a1 || a1 != a2 || epsilon < Matrix<T>::value_epsilon)
    {
        throw std::domain_error("Non-square matrix or epsilon is too small");
    }
}

template <typename T>
static void find_decomposition(const Matrix<T> &a, Matrix<T> &q, Matrix<T> &r)
{
    const typename Matrix<T>::size_type m = a.size1();
    r = a;
    q.identity(m);
    Matrix<T> v, h_k;
    for (typename Matrix<T>::size_type i = 0u; i < m - 1u; ++i)
    {
        v = get_v(r, i);
        h_k = householder_matrix(v);
        q *= h_k;
        r = h_k * r;
    }
}

template <typename T>
static Matrix<T> householder_matrix(const Matrix<T> &v)
{
    const Matrix<T> v_t = transpose(v);
    Matrix<T> result = Matrix<T>().identity(v.size1()) -
        2.0 / (v_t * v)(0u) * v * v_t;
    return result;
}

template <typename T>
static Matrix<T> get_v(const Matrix<T> &a_k,
    const typename Matrix<T>::size_type idx)
{
    const typename Matrix<T>::size_type m = a_k.size1();
    Matrix<T> v(m, 1u);

    typename Matrix<T>::floating_point sum = 0.0;
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        if (i < idx)
        {
            v(i) = 0.0;
        }
        else
        {
            v(i) = a_k(i, idx);
            sum += a_k(i, idx) * a_k(i, idx);
        }
    }
    v(idx) += std::copysign(std::sqrt(sum), a_k(idx, idx));

    return v;
}

template <typename T>
std::pair<
    std::complex<typename Matrix<T>::floating_point>,
    std::complex<typename Matrix<T>::floating_point>>
quadratic_equation(
    const typename Matrix<T>::floating_point a,
    const typename Matrix<T>::floating_point b,
    const typename Matrix<T>::floating_point c)
{
    std::complex<typename Matrix<T>::floating_point> d = b * b - 4.0 * a * c;
    std::complex<typename Matrix<T>::floating_point> x1, x2;
    x1 = (-b + std::sqrt(d)) / a;
    x1 /= 2.0;
    x2 = (-b - std::sqrt(d)) / a;
    x2 /= 2.0;

    return make_pair(x1, x2);
}

template <typename T>
static bool criteria(const Matrix<T> &a_k, const typename Matrix<T>::floating_point epsilon)
{
    const typename Matrix<T>::size_type m = a_k.size1();
    typename Matrix<T>::floating_point sum = 0.0;
    for (typename Matrix<T>::size_type i = 0u; i < m - 1u; ++i)
    {
        for (typename Matrix<T>::size_type j = i + 2u; j < m; ++j)
        {
            sum += a_k(j, i) * a_k(j, i);
        }
    }

    return std::sqrt(sum) > epsilon;
}

#endif
