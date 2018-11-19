#ifndef __HYPERBOLIC_PDE_HPP__
#define __HYPERBOLIC_PDE_HPP__

#include <cstddef>

#include <array>
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

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void to_three_diagonal(ublas::vector<T> &, ublas::vector<T> &, ublas::vector<T> &,
    ublas::vector<T> &, const T, const T) noexcept;

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
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau * tau / (h * h),
        p = 1.0 - b * h / (2.0 * a * a), q = 1.0 + b * h / (2.0 * a * a),
        b_0 = beta - alpha / (h * p) - alpha * h / (2.0 * a * a * tau * tau * p) -
            alpha * d * h / (2.0 * a * a * tau * p) +
            alpha * c * h / (2.0 * a * a * p),
        c_0 = alpha / (h * p),
        a_n = -gamma / (h * q),
        b_n = delta + gamma / (h * q) + gamma * h / (2.0 * a * a * tau * tau * q) +
            gamma * d * h / (2.0 * a * a * tau * q) -
            gamma * c * h / (2.0 * a * a * q);
    if (n_upper < 2U || k_upper < 2U || h < EPSILON || tau < EPSILON || sigma >= 1.0 ||
        !is_enum_includes(initial) || !is_enum_includes(boundary))
    {
        throw std::logic_error("n < 2 || k < 2 || h < epsilon || tau < epsilon || sigma >= 1.0");
    }

    std::array<ublas::vector<T>, THREE> w_h_tau;
    w_h_tau.fill(ublas::vector<T>(n_upper + 1));
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
                tau * tau * f_x_t(j * h, k * tau);
            u_k[j] /= (1.0 + d * tau / 2.0);
        }
        switch (boundary)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                u_k[0] = -alpha / (beta * h - alpha) * u_k[1]
                    + phi_0_t(a, b, c, k * tau) * h / (beta * h - alpha);
                u_k[n_upper] = gamma / (delta * h + gamma) * u_k[n_upper - 1]
                    + phi_l_t(a, b, c, k * tau) * h / (delta * h + gamma);
                break;
            case TWO_POINT_SECOND_ORDER:
                u_k[0] = -c_0 * u_k[1] / b_0 -
                    (alpha * h / (a * a * tau * tau * p) +
                        alpha * d * h / (2.0 * a * a * tau * p)) / b_0 * u_k_minus_1[0] +
                    alpha * h / (2.0 * a * a * tau * tau * p) / b_0 * u_k_minus_2[0] +
                    phi_0_t(a, b, c, k * tau) / b_0 -
                    alpha * h / (2.0 * a * a * p) / b_0 * f_x_t(0.0, k * tau);
                u_k[n_upper] = -a_n * u_k[n_upper - 1] / b_n +
                    (gamma * h / (a * a * tau * tau * q) +
                        gamma * d * h / (2.0 * a * a * tau * q)) / b_n * u_k_minus_1[n_upper] -
                    gamma * h / (2.0 * a * a * tau * tau * q) / b_n * u_k_minus_2[n_upper] +
                    phi_l_t(a, b, c, k * tau) / b_n +
                    gamma * h / (2.0 * a * a * q) / b_n * f_x_t(l, k * tau);
                break;
            case THREE_POINT_SECOND_ORDER:
                u_k[0] =
                    (-4.0 * alpha * u_k[1] + alpha * u_k[2U] +
                        2.0 * phi_0_t(a, b, c, k * tau) * h)
                    / (2.0 * beta * h - 3.0 * alpha);
                u_k[n_upper] =
                    (4.0 * gamma * u_k[n_upper - 1] - gamma * u_k[n_upper - 2U] +
                        2.0 * phi_l_t(a, b, c, k * tau) * h) /
                    (2.0 * delta * h + 3.0 * gamma);
                break;
        }
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
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau * tau / (h * h),
        p = 1.0 - b * h / (2.0 * a * a), q = 1.0 + b * h / (2.0 * a * a);
    if (n_upper < 4U || k_upper < 2U || h < EPSILON || tau < EPSILON ||
        !is_enum_includes(initial) || !is_enum_includes(boundary))
    {
        throw std::logic_error("n < 4 || k < 2 || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper + 1),
        a_j(n_upper + 1, sigma - b * tau * tau / (2.0 * h)),
        b_j(n_upper + 1, c * tau * tau - 2.0 * sigma - 1.0 - d * tau / 2.0),
        c_j(n_upper + 1, sigma + b * tau * tau / (2.0 * h));
    a_j[0] = c_j[n_upper] = 0.0;
    switch (boundary)
    {
        default:
        case TWO_POINT_FIRST_ORDER:
            b_j[0] = beta * h - alpha;
            c_j[0] = alpha;
            a_j[n_upper] = -gamma;
            b_j[n_upper] = gamma + delta * h;
            break;
        case TWO_POINT_SECOND_ORDER:
            b_j[0] = beta - alpha / (h * p) - alpha * h / (2.0 * a * a * tau * tau * p) -
                alpha * d * h / (2.0 * a * a * tau * p) +
                alpha * c * h / (2.0 * a * a * p);
            c_j[0] = alpha / (h * p);
            a_j[n_upper] = -gamma / (h * q);
            b_j[n_upper] = delta + gamma / (h * q) + gamma * h / (2.0 * a * a * tau * tau * q) +
                gamma * d * h / (2.0 * a * a * tau * q) -
                gamma * c * h / (2.0 * a * a * q);
            break;
        case THREE_POINT_SECOND_ORDER:
            break;
    }
    std::array<ublas::vector<T>, THREE> w_h_tau;
    w_h_tau.fill(ublas::vector<T>(n_upper + 1));
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
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            d_j[j] = -2.0 * u_k_minus_1[j] +
                (1.0 - d * tau / 2.0) * u_k_minus_2[j] -
                tau * tau * f_x_t(j * h, k * tau);
        }
        switch (boundary)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                d_j[0] = phi_0_t(a, b, c, k * tau) * h;
                d_j[n_upper] = phi_l_t(a, b, c, k * tau) * h;
                break;
            case TWO_POINT_SECOND_ORDER:
                d_j[0] = -(alpha * h / (a * a * tau * tau * p) +
                        alpha * d * h / (2.0 * a * a * tau * p)) * u_k_minus_1[0] +
                    alpha * h / (2.0 * a * a * tau * tau * p) * u_k_minus_2[0] +
                    phi_0_t(a, b, c, k * tau) -
                    alpha * h / (2.0 * a * a * p) * f_x_t(0.0, k * tau);
                d_j[n_upper] = (gamma * h / (a * a * tau * tau * q) +
                        gamma * d * h / (2.0 * a * a * tau * q)) * u_k_minus_1[n_upper] -
                    gamma * h / (2.0 * a * a * tau * tau * q) * u_k_minus_2[n_upper] +
                    phi_l_t(a, b, c, k * tau) +
                    gamma * h / (2.0 * a * a * q) * f_x_t(l, k * tau);
                break;
            case THREE_POINT_SECOND_ORDER:
                b_j[0] = 2.0 * beta * h - 3.0 * alpha;
                c_j[0] = 4.0 * alpha;
                d_j[0] = 2.0 * phi_0_t(a, b, c, k * tau) * h;
                a_j[n_upper] = -4.0 * gamma;
                b_j[n_upper] = 2.0 * delta * h + 3.0 * gamma;
                d_j[n_upper] = 2.0 * phi_l_t(a, b, c, k * tau) * h;
                to_three_diagonal(a_j, b_j, c_j, d_j, -alpha, gamma);
                break;
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

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void to_three_diagonal(ublas::vector<T> &a, ublas::vector<T> &b, ublas::vector<T> &c,
    ublas::vector<T> &d, const T m_0_2, const T m_n_n_minus_2) noexcept
{
    const std::size_t n = d.size() - 1;
    const T alpha_0 = -m_0_2 / c[1], alpha_n = -m_n_n_minus_2 / a[n - 1];
    b[0] = std::fma(alpha_0, a[1], b[0]);
    c[0] = std::fma(alpha_0, b[1], c[0]);
    d[0] = std::fma(alpha_0, d[1], d[0]);
    a[n] = std::fma(alpha_n, b[n - 1], a[n]);
    b[n] = std::fma(alpha_n, c[n - 1], b[n]);
    d[n] = std::fma(alpha_n, d[n - 1], d[n]);
}

#endif
