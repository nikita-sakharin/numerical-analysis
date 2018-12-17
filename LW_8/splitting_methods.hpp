#ifndef __SPLITTING_METHODS_HPP__
#define __SPLITTING_METHODS_HPP__

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

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::matrix<T> alternating_direction_method(const T, const T, const T,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T, const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<void (const ublas::matrix<T> &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::matrix<T> fractional_step_method(const T, const T, const T,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T, const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<void (const ublas::matrix<T> &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void initial_approximation(const T, const T, const T, std::size_t, std::size_t,
    const T, const T, const std::function<T (const T &, const T &, const T &,
        const T &, const T &)> &, ublas::matrix<T> &);

template<typename T, typename>
ublas::matrix<T> alternating_direction_method(const T a, const T b, const T mu,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &f_x_y_t,
    const T l_1, const std::size_t n_1, const T l_2, const std::size_t n_2,
    const T t, const std::size_t k_upper,
    const T alpha_1, const T beta_1, const T alpha_2, const T beta_2,
    const T alpha_3, const T beta_3, const T alpha_4, const T beta_4,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_1_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_2_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_3_y_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_4_y_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &psi_x_y,
    const std::function<void (const ublas::matrix<T> &)> &get_error)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h_1 = l_1 / n_1, h_2 = l_2 / n_2, tau = t / k_upper;
    if (n_1 < 2U || n_2 < 2U || !k_upper ||
        h_1 < EPSILON || h_1 < EPSILON, tau / 2.0 < EPSILON)
    {
        throw std::logic_error("n_1 < 2 || n_2 < 2 || !k || "
            "h_1 < epsilon || h_2 < epsilon || tau < epsilon");
    }

    ublas::matrix<T> u_k(n_1 + 1, n_2 + 1), u_k_minus_1_divides_2(n_1 + 1, n_2 + 1),
        u_k_minus_1(n_1 + 1, n_2 + 1);
    initial_approximation(a, b, mu, n_1, n_2, h_1, h_2, psi_x_y, u_k_minus_1);
    get_error(u_k_minus_1);
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            u_k[j] = (sigma + b * tau / (2.0 * h)) * u_k_minus_1[j + 1] +
                (1.0 - 2.0 * sigma + c * tau) * u_k_minus_1[j] +
                (sigma - b * tau / (2.0 * h)) * u_k_minus_1[j - 1];
        }
        u_k[0] =
            (-alpha * u_k[1] + phi_0_t(a, b, c, k * tau) * h) /
            (beta * h - alpha);
        u_k[n_upper] =
            (gamma * u_k[n_upper - 1] + phi_l_t(a, b, c, k * tau) * h) /
            (delta * h + gamma);
        u_k.swap(u_k_minus_1);
        get_error(u_k);
    }

    return u_k_minus_1;
}

template<typename T, typename>
ublas::matrix<T> fractional_step_method(const T a, const T b, const T mu,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &f_x_y_t,
    const T l_1, const std::size_t n_1, const T l_2, const std::size_t n_2,
    const T t, const std::size_t k,
    const T alpha_1, const T beta_1, const T alpha_2, const T beta_2,
    const T alpha_3, const T beta_3, const T alpha_4, const T beta_4,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_1_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_2_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_3_y_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_4_y_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &psi_x_y,
    const std::function<void (const ublas::matrix<T> &)> &get_error)
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
static void initial_approximation(const T a, const T b, const T mu,
    const std::size_t n_1, const std::size_t n_2, const T h_1, const T h_2,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &)> &psi_x_y, ublas::matrix<T> &u_0)
{
    for (std::size_t i = 0; i <= n_1; ++i)
    {
        for (std::size_t j = 0; j <= n_2; ++j)
        {
            u_0(i, j) = psi_x_y(a, b, mu, i * h_1, j * h_2);
        }
    }
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
