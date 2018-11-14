#ifndef __HYPERBOLIC_PDE_HPP__
#define __HYPERBOLIC_PDE_HPP__

#include <cstddef>

#include <limits>
#include <type_traits>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../generic/thomas_algorithm.hpp"

enum NumDiff : uint
{
    TWO_POINT_FIRST_ORDER,
    TWO_POINT_SECOND_ORDER,
    THREE_POINT_SECOND_ORDER,
};

static constexpr std::size_t THREE = 3U;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> explicit_fdm(const T, const T, const T, const T,
    const std::function<T (const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    NumDiff, NumDiff);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> implicit_fdm(const T, const T, const T, const T,
    const std::function<T (const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    NumDiff, NumDiff);

static bool is_enum_includes(NumDiff) noexcept;

template<typename T, typename>
ublas::vector<T> explicit_fdm(const T a, const T b, const T c, const T d,
    const std::function<T (const T &, const T &)> &f_x_t,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_1_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_prime_1_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_prime_prime_1_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_2_x,
    const NumDiff initial, const NumDiff boundary)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau * tau / (h * h);
    if (!n_upper || !k_upper || h < EPSILON || tau < EPSILON || sigma >= 1.0 ||
        !is_enum_includes(initial) || !is_enum_includes(boundary))
    {
        throw std::logic_error("!n || !k || h < epsilon || tau < epsilon || sigma >= 1.0");
    }

    ublas::vector<ublas::vector<T>> w_h_tau(THREE, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_1_x(a, b, c, j * h);
        switch (initial)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                w_h_tau[1][j] = w_h_tau[0][j] + psi_2_x(a, b, c, j * h) * tau;
                break;
            case TWO_POINT_SECOND_ORDER:
                w_h_tau[1][j] = (1.0 + c * tau * tau / 2.0) * w_h_tau[0][j] +
                    (tau - d * tau * tau / 2.0) * psi_2_x(a, b, c, j * h) +
                    a * a * tau * tau / 2.0 * psi_prime_prime_1_x(a, b, c, j * h) +
                    b * tau * tau / 2.0 * psi_prime_1_x(a, b, c, j * h) +
                    tau * tau / 2.0 * f_x_t(j * h, 0.0);
                break;
            case THREE_POINT_SECOND_ORDER:
                throw std::logic_error("initial: THREE_POINT_SECOND_ORDER");
        }
    }
    for (std::size_t k = 2U; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % THREE],
            &u_k_minus_1 = w_h_tau[(k - 1) % THREE],
            &u_k_minus_2 = w_h_tau[(k - 2U) % THREE];
        for (std::size_t j = 1U; j < n_upper; ++j)
        {
            u_k[j] = (sigma + b * tau * tau / (2.0 * h)) * u_k_minus_1[j + 1] +
                (2.0 - 2.0 * sigma + c * tau * tau) * u_k_minus_1[j] +
                (sigma - b * tau * tau / (2.0 * h)) * u_k_minus_1[j - 1] +
                (d * tau / 2.0 - 1.0) * u_k_minus_2[j] +
                tau * tau * f_x_t(j * h, tau * (k + 1));
            u_k[j] /= (1.0 + d * tau / 2.0);
        }
        u_k[0] = -alpha / (beta * h - alpha) * u_k[1]
            + phi_0_t(a, b, c, k * tau) * h / (beta * h - alpha);
        u_k[n_upper] = gamma / (delta * h + gamma) * u_k[n_upper - 1]
            + phi_l_t(a, b, c, k * tau) * h / (delta * h + gamma);
    }

    return w_h_tau[k_upper % THREE];
}

template<typename T, typename>
ublas::vector<T> implicit_fdm(const T a, const T b, const T c, const T d,
    const std::function<T (const T &, const T &)> &f_x_t,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_1_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_prime_1_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_prime_prime_1_x,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_2_x,
    const NumDiff initial, const NumDiff boundary)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau * tau / (h * h);
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON ||
        !is_enum_includes(initial) || !is_enum_includes(boundary))
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper + 1),
        a_j(n_upper + 1, sigma - b * tau * tau / (2.0 * h)),
        b_j(n_upper + 1, c * tau * tau - 2.0 * sigma - 1.0 - d * tau / 2.0),
        c_j(n_upper + 1, sigma + b * tau * tau / (2.0 * h));
    a_j[0] = c_j[n_upper] = 0.0;
    b_j[0] = beta * h - alpha;
    c_j[0] = alpha;
    a_j[n_upper] = -gamma;
    b_j[n_upper] = gamma + delta * h;
    ublas::vector<ublas::vector<T>> w_h_tau(THREE, ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_1_x(a, b, c, j * h);
        switch (initial)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                w_h_tau[1][j] = w_h_tau[0][j] + psi_2_x(a, b, c, j * h) * tau;
                break;
            case TWO_POINT_SECOND_ORDER:
                w_h_tau[1][j] = (1.0 + c * tau * tau / 2.0) * w_h_tau[0][j] +
                    (tau - d * tau * tau / 2.0) * psi_2_x(a, b, c, j * h) +
                    a * a * tau * tau / 2.0 * psi_prime_prime_1_x(a, b, c, j * h) +
                    b * tau * tau / 2.0 * psi_prime_1_x(a, b, c, j * h) +
                    tau * tau / 2.0 * f_x_t(j * h, 0.0);
                break;
            case THREE_POINT_SECOND_ORDER:
                throw std::logic_error("initial: THREE_POINT_SECOND_ORDER");
        }
    }
    for (std::size_t k = 2U; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % THREE],
            &u_k_minus_1 = w_h_tau[(k - 1) % THREE],
            &u_k_minus_2 = w_h_tau[(k - 2U) % THREE];
        d_j[0] = phi_0_t(a, b, c, k * tau) * h;
        d_j[n_upper] = phi_l_t(a, b, c, k * tau) * h;
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            d_j[j] = -2.0 * u_k_minus_1[j] +
                (1.0 - d * tau / 2.0) * u_k_minus_2[j] -
                tau * tau * f_x_t(j * h, k * tau);
        }
        u_k = thomas_algorithm(a_j, b_j, c_j, d_j);
    }

    return w_h_tau[k_upper % THREE];
}

static bool is_enum_includes(const NumDiff num_diff) noexcept
{
    switch (num_diff)
    {
        case TWO_POINT_FIRST_ORDER:
        case TWO_POINT_SECOND_ORDER:
        case THREE_POINT_SECOND_ORDER:
            return true;
        default:
            return false;
    }
}

#endif
