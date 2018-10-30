#ifndef __PARABOLIC_PDE_HPP__
#define __PARABOLIC_PDE_HPP__

#include <cstddef>

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
ublas::vector<T> explicit_fdm(const T, const T, std::size_t,
    const T, std::size_t, const std::function<T (const T &, const T &)> &,
    const std::function<T (const T &, const T &)> &,
    const std::function<T (const T &, const T &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> implicit_fdm(const T, const T, std::size_t,
    const T, std::size_t, const std::function<T (const T &, const T &)> &,
    const std::function<T (const T &, const T &)> &,
    const std::function<T (const T &, const T &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> crank_nicolson(const T, const T, std::size_t,
    const T, std::size_t, const std::function<T (const T &, const T &)> &,
    const std::function<T (const T &, const T &)> &,
    const std::function<T (const T &, const T &)> &);

template<typename T, typename>
ublas::vector<T> explicit_fdm(const T a, const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const std::function<T (const T &, const T &)> &u_0_t,
    const std::function<T (const T &, const T &)> &u_l_t,
    const std::function<T (const T &, const T &)> &u_x_0)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    if (!n_upper || !k_upper || h < EPSILON || tau < EPSILON || sigma > 0.5)
    {
        throw std::logic_error("!n || !k || h < epsilon || tau < epsilon || sigma > 0.5");
    }

    ublas::vector<ublas::vector<T>> w_h_tau(TWO, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = u_x_0(a, j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        u_k[0]       = u_0_t(a, k * tau);
        u_k[n_upper] = u_l_t(a, k * tau);
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            u_k[j] = sigma * u_k_minus_1[j + 1] +
                (1.0 - 2.0 * sigma) * u_k_minus_1[j] +
                sigma * u_k_minus_1[j - 1];
        }
    }

    return w_h_tau[k_upper % TWO];
}

template<typename T, typename>
ublas::vector<T> implicit_fdm(const T a, const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const std::function<T (const T &, const T &)> &u_0_t,
    const std::function<T (const T &, const T &)> &u_l_t,
    const std::function<T (const T &, const T &)> &u_x_0)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    const std::size_t n_upper_minus_1 = n_upper - 1;
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper_minus_1),
        a_j(n_upper_minus_1, sigma),
        b_j(n_upper_minus_1, -(1.0 + 2.0 * sigma)),
        c_j(n_upper_minus_1, sigma);
    a_j[0] = c_j[n_upper_minus_1 - 1] = 0.0;
    ublas::vector<ublas::vector<T>> w_h_tau(TWO, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = u_x_0(a, j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        u_k[0]       = u_0_t(a, k * tau);
        u_k[n_upper] = u_l_t(a, k * tau);
        d_j[0] =
            -(u_k_minus_1[1] + sigma * u_0_t(a, k * tau));
        d_j[n_upper_minus_1 - 1] =
            -(u_k_minus_1[n_upper - 1] + sigma * u_l_t(a, k * tau));
        for (std::size_t j = 1; j < n_upper_minus_1 - 1; ++j)
        {
            d_j[j] = -u_k_minus_1[j + 1];
        }
        x = thomas_algorithm(a_j, b_j, c_j, d_j);
        std::move(x.cbegin(), x.cend(), u_k.begin() + 1);
    }

    return w_h_tau[k_upper % TWO];
}

template<typename T, typename>
ublas::vector<T> crank_nicolson(const T a, const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const std::function<T (const T &, const T &)> &u_0_t,
    const std::function<T (const T &, const T &)> &u_l_t,
    const std::function<T (const T &, const T &)> &u_x_0)
{
    static constexpr T THETA = 0.5, EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    const std::size_t n_upper_minus_1 = n_upper - 1;
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper_minus_1),
        a_j(n_upper_minus_1, THETA * sigma),
        b_j(n_upper_minus_1, -(1.0 + THETA * 2.0 * sigma)),
        c_j(n_upper_minus_1, THETA * sigma);
    a_j[0] = c_j[n_upper_minus_1 - 1] = 0.0;
    ublas::vector<ublas::vector<T>> w_h_tau(TWO, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = u_x_0(a, j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        u_k[0]       = u_0_t(a, k * tau);
        u_k[n_upper] = u_l_t(a, k * tau);
        d_j[0] =
            -(u_k_minus_1[1] + THETA * sigma * u_0_t(a, k * tau)) -
            (1.0 - THETA) * sigma *
            (u_k_minus_1[2U] - 2.0 * u_k_minus_1[1] + u_k_minus_1[0]);
        d_j[n_upper_minus_1 - 1] =
            -(u_k_minus_1[n_upper - 1] + THETA * sigma * u_l_t(a, k * tau)) -
            (1.0 - THETA) * sigma *
            (u_k_minus_1[n_upper - 2U] - 2.0 * u_k_minus_1[n_upper - 1] + u_k_minus_1[n_upper]);
        for (std::size_t j = 1; j < n_upper_minus_1 - 1; ++j)
        {
            d_j[j] = (THETA - 1.0) * sigma *
                (u_k_minus_1[j + 2U] - 2.0 * u_k_minus_1[j + 1] + u_k_minus_1[j]) -
                u_k_minus_1[j + 1];
        }
        x = thomas_algorithm(a_j, b_j, c_j, d_j);
        std::move(x.cbegin(), x.cend(), u_k.begin() + 1);
    }

    return w_h_tau[k_upper % TWO];
}

#endif
