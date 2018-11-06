#ifndef __PARABOLIC_PDE_HPP__
#define __PARABOLIC_PDE_HPP__

#include <cstddef>

#include <limits>
#include <type_traits>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../generic/thomas_algorithm.hpp"

static constexpr std::size_t TWO = 2U;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> explicit_fdm(const T, const T, const T,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> implicit_fdm(const T, const T, const T,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> crank_nicolson(const T, const T, const T,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &);

template<typename T, typename>
ublas::vector<T> explicit_fdm(const T a, const T b, const T c,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_x)
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
        w_h_tau[0][j] = psi_x(a, b, c, j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            u_k[j] = sigma * u_k_minus_1[j + 1] +
                (1.0 - 2.0 * sigma) * u_k_minus_1[j] +
                sigma * u_k_minus_1[j - 1] +
                b * tau * (u_k_minus_1[j + 1] - u_k_minus_1[j - 1]) / (2.0 * h) +
                c * tau * u_k_minus_1[j];
        }
        u_k[0] = -alpha / (beta * h - alpha) * u_k[1]
            + phi_0_t(a, b, c, k * tau) * h / (beta * h - alpha);
        u_k[n_upper] = gamma / (delta * h + gamma) * u_k[n_upper - 1]
            + phi_l_t(a, b, c, k * tau) * h / (delta * h + gamma);
    }

    return w_h_tau[k_upper % TWO];
}

template<typename T, typename>
ublas::vector<T> implicit_fdm(const T a, const T b, const T c,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_x)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper + 1),
        a_j(n_upper + 1, sigma - b * tau / (2.0 * h)),
        b_j(n_upper + 1, c * tau - (1.0 + 2.0 * sigma)),
        c_j(n_upper + 1, sigma + b * tau / (2.0 * h));
    a_j[0] = c_j[n_upper] = 0.0;
    b_j[0] = beta * h - alpha;
    c_j[0] = alpha;
    a_j[n_upper] = -gamma;
    b_j[n_upper] = gamma + delta * h;
    ublas::vector<ublas::vector<T>> w_h_tau(TWO, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_x(a, b, c, j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        d_j[0] = phi_0_t(a, b, c, k * tau) * h;
        d_j[n_upper] = phi_l_t(a, b, c, k * tau) * h;
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            d_j[j] = -u_k_minus_1[j];
        }
        u_k = thomas_algorithm(a_j, b_j, c_j, d_j);
    }

    return w_h_tau[k_upper % TWO];
}

template<typename T, typename>
ublas::vector<T> crank_nicolson(const T a, const T b, const T c,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_x)
{
    static constexpr T THETA = 0.5, EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper + 1),
        a_j(n_upper + 1, THETA * (sigma - b * tau / (2.0 * h))),
        b_j(n_upper + 1, -1.0 + THETA * (c * tau - 2.0 * sigma)),
        c_j(n_upper + 1, THETA * (sigma + b * tau / (2.0 * h)));
    a_j[0] = c_j[n_upper] = 0.0;
    b_j[0] = beta * h - alpha;
    c_j[0] = alpha;
    a_j[n_upper] = -gamma;
    b_j[n_upper] = gamma + delta * h;
    ublas::vector<ublas::vector<T>> w_h_tau(TWO, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_x(a, b, c, j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        d_j[0] = phi_0_t(a, b, c, k * tau) * h;
        d_j[n_upper] = phi_l_t(a, b, c, k * tau) * h;
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            d_j[j] = -u_k_minus_1[j] + (THETA - 1.0) *
                (sigma * (u_k_minus_1[j + 1] - 2.0 * u_k_minus_1[j] + u_k_minus_1[j - 1]) +
                b * tau / (2.0 * h) * (u_k_minus_1[j + 1] - u_k_minus_1[j - 1]) +
                c * tau * u_k_minus_1[j]);
        }
        u_k = thomas_algorithm(a_j, b_j, c_j, d_j);
    }

    return w_h_tau[k_upper % TWO];
}

#endif
