#ifndef __JACOBI_EIGENVALUE_ALGORITHM_HPP__
#define __JACOBI_EIGENVALUE_ALGORITHM_HPP__

#include <utility>

#include "./matrix/matrix.hpp"

template <typename T>
std::pair<Matrix<T>, Matrix<T>> jacobi_eigenvalue_algorithm(const Matrix<T> &,
    const typename Matrix<T>::floating_point);

template <typename T>
static inline void arg_vilid(typename Matrix<T>::size_type, typename Matrix<T>::size_type,
    typename Matrix<T>::floating_point);

template <typename T>
static std::pair<typename Matrix<T>::size_type, typename Matrix<T>::size_type>
    find_max(const Matrix<T> &);

template <typename T>
static typename Matrix<T>::floating_point t(const Matrix<T> &);

template <typename T>
static Matrix<T> rotation_matrix(const Matrix<T> &, const std::pair<typename Matrix<T>::size_type,
    typename Matrix<T>::size_type> &);

template <typename T>
std::pair<Matrix<T>, Matrix<T>> jacobi_eigenvalue_algorithm(const Matrix<T> &a,
    const typename Matrix<T>::floating_point epsilon)
{
    arg_vilid<T>(a.size1(), a.size2(), epsilon);
    Matrix<T> a_k = a, u_k, mult;
    mult.identity(a.size1());
    while (t(a_k) > epsilon)
    {
        std::pair<typename Matrix<T>::size_type,
            typename Matrix<T>::size_type> max = find_max(a_k);
        u_k = rotation_matrix(a_k, max);
        a_k = transpose(u_k) * a_k * u_k;
        mult *= u_k;
    }

    return std::make_pair(a_k, mult);
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
static std::pair<typename Matrix<T>::size_type, typename Matrix<T>::size_type>
    find_max(const Matrix<T> &a)
{
    typename Matrix<T>::size_type i_max, j_max;
    i_max = j_max = Matrix<T>::size_max;
    const typename Matrix<T>::size_type m = a.size1(), n = a.size2();
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        for (typename Matrix<T>::size_type j = i + 1u; j < n; ++j)
        {
            if (i_max == Matrix<T>::size_max || std::fabs(a(i_max, j_max)) < std::fabs(a(i, j)))
            {
                i_max = i;
                j_max = j;
            }
        }
    }

    return std::make_pair(i_max, j_max);
}

template <typename T>
static Matrix<T> rotation_matrix(const Matrix<T> &a,
    const std::pair<typename Matrix<T>::size_type, typename Matrix<T>::size_type> &max)
{
    const typename Matrix<T>::size_type i = max.first, j = max.second;
    typename Matrix<T>::floating_point
        phi_k = 0.5 * std::atan(2.0 * a(i, j) / (a(i, i) - a(j, j)));
    Matrix<T> result;
    result.identity(a.size1());
    result(i, i) = result(j, j) = std::cos(phi_k);
    result(i, j) = -std::sin(phi_k);
    result(j, i) = std::sin(phi_k);

    return result;
}

template <typename T>
static typename Matrix<T>::floating_point t(const Matrix<T> &a)
{
    const typename Matrix<T>::size_type m = a.size1(), n = a.size2();
    typename Matrix<T>::floating_point sum = 0.0;
    for (typename Matrix<T>::size_type i = 0u; i < m; ++i)
    {
        for (typename Matrix<T>::size_type j = i + 1u; j < n; ++j)
        {
            sum += a(i, j) * a(i, j);
        }
    }

    return std::sqrt(sum);
}

#endif
