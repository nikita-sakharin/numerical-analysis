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
ublas::matrix<T> successive_fdm(const T, const T, const T,
    const std::function<T (const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T, const T, const T, const T, const T,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    const std::function<void (const ublas::matrix<T> &)> &, const T);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::matrix<T> seidel_fdm(const T, const T, const T,
    const std::function<T (const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T, const T, const T, const T, const T,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    const T, const std::function<void (const ublas::matrix<T> &)> &, const T);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void initial_approximation(std::size_t, std::size_t, std::size_t, std::size_t,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    ublas::matrix<T> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static T norm_assign(T &, const T, const T) noexcept;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static T next_u_k(const T, const T, const T,
    const std::function<T (const T &, const T &)> &, const T, const T,
    const T, const T, const T, const T, const T, const T, const T, const T,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &,
    const T, ublas::matrix<T> &, const ublas::matrix<T> &) noexcept;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static constexpr T relax(const T omega,
    const T u_k_i_j, const T u_k_minus_1_i_j) noexcept;

template<typename T, typename>
ublas::matrix<T> successive_fdm(const T a, const T b, const T c,
    const std::function<T (const T &, const T &)> &f_x_y,
    const T l_1, const std::size_t n_1,
    const T l_2, const std::size_t n_2,
    const T alpha_1, const T beta_1, const T alpha_2, const T beta_2,
    const T alpha_3, const T beta_3, const T alpha_4, const T beta_4,
    const std::function<T (const T &)> &phi_1_y,
    const std::function<T (const T &)> &phi_2_y,
    const std::function<T (const T &)> &phi_3_x,
    const std::function<T (const T &)> &phi_4_x,
    const std::function<void (const ublas::matrix<T> &)> &get_error,
    const T epsilon)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h_1 = l_1 / n_1, h_2 = l_2 / n_2;
    if (n_1 < 2U || n_2 < 2U || h_1 < EPSILON || h_1 < EPSILON)
    {
        throw std::logic_error("n_1 < 2 || n_2 < 2 || h_1 < epsilon || h_2 < epsilon");
    }

    ublas::matrix<T> u_k(n_1 + 1, n_2 + 1), u_k_minus_1(n_1 + 1, n_2 + 1);
    initial_approximation(n_1, n_2, h_1, h_2, phi_1_y, phi_2_y, u_k_minus_1);
    get_error(u_k_minus_1);
    T norm;
    do
    {
        norm = next_u_k<T>(a, b, c, f_x_y, h_1, h_2, alpha_1, beta_1, alpha_2, beta_2,
            alpha_3, beta_3, alpha_4, beta_4,
            phi_1_y, phi_2_y, phi_3_x, phi_4_x, 1.0, u_k, u_k_minus_1);
        get_error(u_k);
        u_k.swap(u_k_minus_1);
    } while (norm >= epsilon);

    return u_k_minus_1;
}

template<typename T, typename>
ublas::matrix<T> seidel_fdm(const T a, const T b, const T c,
    const std::function<T (const T &, const T &)> &f_x_y,
    const T l_1, const std::size_t n_1,
    const T l_2, const std::size_t n_2,
    const T alpha_1, const T beta_1, const T alpha_2, const T beta_2,
    const T alpha_3, const T beta_3, const T alpha_4, const T beta_4,
    const std::function<T (const T &)> &phi_1_y,
    const std::function<T (const T &)> &phi_2_y,
    const std::function<T (const T &)> &phi_3_x,
    const std::function<T (const T &)> &phi_4_x,
    const T omega, const std::function<void (const ublas::matrix<T> &)> &get_error,
    const T epsilon)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h_1 = l_1 / n_1, h_2 = l_2 / n_2;
    if (n_1 < 2U || n_2 < 2U || h_1 < EPSILON || h_1 < EPSILON || omega <= 0)
    {
        throw std::logic_error("n_1 < 2 || n_2 < 2 || h_1 < epsilon || h_2 < epsilon || omega <= 0");
    }

    ublas::matrix<T> u_k(n_1 + 1, n_2 + 1);
    initial_approximation(n_1, n_2, h_1, h_2, phi_1_y, phi_2_y, u_k);
    get_error(u_k);
    T norm;
    do
    {
        norm = next_u_k<T>(a, b, c, f_x_y, h_1, h_2, alpha_1, beta_1, alpha_2, beta_2,
            alpha_3, beta_3, alpha_4, beta_4,
            phi_1_y, phi_2_y, phi_3_x, phi_4_x, omega, u_k, u_k);
        get_error(u_k);
    } while (norm >= epsilon);

    return u_k;
}

template<typename T, typename>
static void initial_approximation(const std::size_t n_1, const std::size_t n_2,
    const std::size_t h_1, const std::size_t h_2,
    const std::function<T (const T &)> &phi_1_y,
    const std::function<T (const T &)> &phi_2_y,
    ublas::matrix<T> &u_0)
{
    for (std::size_t j = 0; j <= n_2; ++j)
    {
        const T u_0_j = phi_1_y(j * h_2), u_n_1_j = phi_2_y(j * h_2);
        for (std::size_t i = 0; i <= n_1; ++i)
        {
            u_0(i, j) = u_exact(i * h_1, j * h_2);
//          u_0(i, j) = u_0_j + (u_n_1_j - u_0_j) * i / n_1;
            static_cast<void>(u_0_j);
            static_cast<void>(u_n_1_j);
        }
    }
}

template<typename T, typename>
static T norm_assign(T &u_k_i_j, const T value, const T norm) noexcept
{
    const T new_norm = std::abs(u_k_i_j - value);
    u_k_i_j = value;
    return std::max(norm, new_norm);
}

template<typename T, typename>
static T next_u_k(const T a, const T b, const T c,
    const std::function<T (const T &, const T &)> &f_x_y, const T h_1, const T h_2,
    const T alpha_1, const T beta_1, const T alpha_2, const T beta_2,
    const T alpha_3, const T beta_3, const T alpha_4, const T beta_4,
    const std::function<T (const T &)> &phi_1_y,
    const std::function<T (const T &)> &phi_2_y,
    const std::function<T (const T &)> &phi_3_x,
    const std::function<T (const T &)> &phi_4_x,
    const T omega, ublas::matrix<T> &u_k, const ublas::matrix<T> &u_k_minus_1) noexcept
{
    const T div = (h_1 * h_1 * h_2 * h_2 * c + 2.0 * (h_1 * h_1 + h_2 * h_2));
    const std::size_t n_1 = u_k.size1() - 1, n_2 = u_k.size2() - 1;
    T u_k_i_j, norm = 0.0;
    for (std::size_t i = 1; i < n_1; ++i)
    {
        for (std::size_t j = 1; j < n_2; ++j)
        {
            u_k_i_j = (
                (h_2 * h_2 - a / 2.0 * h_1 * h_2 * h_2) * u_k_minus_1(i + 1, j) +
                (h_2 * h_2 + a / 2.0 * h_1 * h_2 * h_2) * u_k_minus_1(i - 1, j) +
                (h_1 * h_1 - b / 2.0 * h_1 * h_1 * h_2) * u_k_minus_1(i, j + 1) +
                (h_1 * h_1 + b / 2.0 * h_1 * h_1 * h_2) * u_k_minus_1(i, j - 1) -
                h_1 * h_1 * h_2 * h_2 * f_x_y(i * h_1, j * h_2)) / div;
            u_k_i_j = relax(omega, u_k_i_j, u_k(i, j));
            norm = norm_assign(u_k(i, j), u_k_i_j, norm);
        }
    }
    for (std::size_t j = 1; j < n_2; ++j)
    {
        u_k_i_j = phi_1_y(j * h_2) * h_1 / (beta_1 * h_1 - alpha_1) -
                alpha_1 * u_k(1, j) / (beta_1 * h_1 - alpha_1);
        norm = norm_assign(u_k(0, j), u_k_i_j, norm);
        u_k_i_j = phi_2_y(j * h_2) * h_1 / (beta_2 * h_1 + alpha_2) +
                alpha_2 * u_k(n_1 - 1, j) / (beta_2 * h_1 + alpha_2);
        norm = norm_assign(u_k(n_1, j), u_k_i_j, norm);
    }
    for (std::size_t i = 1; i < n_1; ++i)
    {
        u_k_i_j = phi_3_x(i * h_1) * h_2 / (beta_3 * h_2 - alpha_3) -
                alpha_3 * u_k(i, 1) / (beta_3 * h_2 - alpha_3);
        norm = norm_assign(u_k(i, 0), u_k_i_j, norm);
        u_k_i_j = phi_4_x(i * h_1) * h_2 / (beta_4 * h_2 + alpha_4) +
                alpha_4 * u_k(i, n_2 - 1) / (beta_4 * h_2 + alpha_4);
        norm = norm_assign(u_k(i, n_2), u_k_i_j, norm);
    }

    return norm;
}

template<typename T, typename>
inline static constexpr T relax(const T omega,
    const T u_k_i_j, const T u_k_minus_1_i_j) noexcept
{
    return (1.0 - omega) * u_k_minus_1_i_j + omega * u_k_i_j;
}

#endif
