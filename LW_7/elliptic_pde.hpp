#ifndef __ELLIPTIC_PDE_HPP__
#define __ELLIPTIC_PDE_HPP__

#include <cstddef>

#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../generic/thomas_algorithm.hpp"

static constexpr std::size_t TWO = 2U;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::matrix<T> liebmann_fdm(const T, const T, const T,
    const std::function<T (const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::matrix<T> seidel_fdm(const T, const T, const T,
    const std::function<T (const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
T norm_max_diff(const ublas::matrix<T> &, const ublas::matrix<T> &) noexcept;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
void initial_approximation(std::size_t, std::size_t, std::size_t, std::size_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    ublas::matrix<T> &);

template<typename T, typename>
ublas::matrix<T> liebmann_fdm(const T a, const T b, const T c,
    const std::function<T (const T &, const T &)> &f_x_y,
    const T l_1, const std::size_t n_1,
    const T l_2, const std::size_t n_2,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_1_y,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_2_y,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_3_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_4_x)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h_1 = l_1 / n_1, h_2 = l_2 / n_2,
        div = (h_1 * h_1 * h_2 * h_2 * c + 2.0 * (h_1 * h_1 + h_2 * h_2));
    if (n_1 < 2U || n_2 < 2U || h_1 < EPSILON || h_1 < EPSILON)
    {
        throw std::logic_error("n_1 < 2 || n_2 < 2 || h_1 < epsilon || h_2 < epsilon");
    }

    std::array<ublas::matrix<T>, TWO> w_h_1_h_2;
    w_h_1_h_2.fill(ublas::matrix<T>(n_1 + 1, n_2 + 1));
    initial_approximation(n_1, n_2, h_1, h_2, phi_1_y, phi_2_y, w_h_1_h_2[0]);
    do
    {
        ublas::vector<T> &u_k = w_h_1_h_2[1],
            &u_k_minus_1 = w_h_1_h_2[0];
        for (std::size_t i = 0; i <= n_1; ++i)
        {
            for (std::size_t j = 0; j <= n_2; ++j)
            {
                u_k(i, j) =
                    (h_2 * h_2 - a / 2.0 * h_1 * h_2 * h_2) * u_k_minus_1(i + 1, j) +
                    (h_2 * h_2 + a / 2.0 * h_1 * h_2 * h_2) * u_k_minus_1(i - 1, j) +
                    (h_1 * h_1 + b / 2.0 * h_1 * h_1 * h_2) * u_k_minus_1(i, j + 1) +
                    (h_1 * h_1 + b / 2.0 * h_1 * h_1 * h_2) * u_k_minus_1(i, j - 1) -
                    h_1 * h_1 * h_2 * h_2 * f_x_y(i * h_1, j * h_2);
                u_k /= div;
            }
        }
        u_k.swap(u_k_minus_1);
    } while (norm_max_diff(w_h_1_h_2[0], w_h_1_h_2[1]) >= epsilon)

    return w_h_1_h_2[0];
}

template<typename T, typename>
ublas::matrix<T> seidel_fdm(const T a, const T b, const T c,
    const std::function<T (const T &, const T &)> &f_x_y,
    const T l_1, const std::size_t n_1,
    const T l_2, const std::size_t n_2,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_1_y,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_2_y,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_3_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_4_x)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h_1 = l_1 / n_1, h_2 = l_2 / n_2,
        div = (h_1 * h_1 * h_2 * h_2 * c + 2.0 * (h_1 * h_1 + h_2 * h_2));
    if (n_1 < 2U || n_2 < 2U || h_1 < EPSILON || h_1 < EPSILON)
    {
        throw std::logic_error("n_1 < 2 || n_2 < 2 || h_1 < epsilon || h_2 < epsilon");
    }

    ublas::matrix<T> w_h_1_h_2(n_1 + 1, n_2 + 1);
    initial_approximation(n_1, n_2, h_1, h_2, phi_1_y, phi_2_y, w_h_1_h_2);
    do
    {
        ublas::vector<T> &u_k = w_h_1_h_2;
        for (std::size_t i = 0; i <= n_1; ++i)
        {
            for (std::size_t j = 0; j <= n_2; ++j)
            {
                const T u_k_minus_1_i_j = u_k(i, j);
                u_k(i, j) =
                    (h_2 * h_2 - a / 2.0 * h_1 * h_2 * h_2) * u_k(i + 1, j) +
                    (h_2 * h_2 + a / 2.0 * h_1 * h_2 * h_2) * u_k(i - 1, j) +
                    (h_1 * h_1 + b / 2.0 * h_1 * h_1 * h_2) * u_k(i, j + 1) +
                    (h_1 * h_1 + b / 2.0 * h_1 * h_1 * h_2) * u_k(i, j - 1) -
                    h_1 * h_1 * h_2 * h_2 * f_x_y(i * h_1, j * h_2);
                u_k /= div;
                norm = std::max(norm, std::abs(u_k(i, j) - u_k_minus_1_i_j));
            }
        }
    } while (norm >= epsilon)

    return w_h_1_h_2;
}

template<typename T, typename>
T norm_max_diff(const ublas::matrix<T> &matrix1, const ublas::matrix<T> &matrix2) noexcept
{
    T norm = 0.0;
    const std::size_t size1 = matrix1.size1(), size2 = matrix1.size2();
    for (std::size_t i = 0; i < size1; ++i)
    {
        for (std::size_t j = 0; j < size2; ++j)
        {
            norm = std::max(norm, std::abs(matrix1(i, j) - matrix2(i, j)));
        }
    }

    return norm;
}

template<typename T, typename>
void initial_approximation(const std::size_t n_1, const std::size_t n_2,
    const std::size_t h_1, const std::size_t h_2,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_1_y,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_2_y,
    ublas::matrix<T> &u_0)
{
    for (std::size_t j = 0; j <= n_2; ++j)
    {
        const T u_0_j = phi_1_y(j * h_2), u_n_1_j = phi_2_y(j * h_2);
        for (std::size_t i = 0; i <= n_1; ++i)
        {
            u_0(i, j) = u_0_j + (u_n_1_j - u_0_j) * i / n_1;
        }
    }
}

#endif
