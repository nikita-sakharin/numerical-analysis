#ifndef __SPLITTING_METHODS_HPP__
#define __SPLITTING_METHODS_HPP__

#include <cstddef>

#include <algorithm>
#include <limits>
#include <type_traits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../generic/thomas_algorithm.hpp"

enum Method : uint
{
    ALTERNATING_DIRECTION,
    FRACTIONAL_STEP
};

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::matrix<T> splitting_method(const T, const T, const T,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &,
    const T, std::size_t, const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T, const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    Method, const std::function<void (const ublas::matrix<T> &)> &);

static constexpr bool is_enum_includes(Method) noexcept;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void initial_approximation(const T, const T, const T,
    std::size_t, const T, std::size_t, const T,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &)> &, ublas::matrix<T> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void step_x(const T, const T, const T,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &, std::size_t, const T, std::size_t, const T,
    const T, std::size_t, const T, const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const ublas::vector<T> &, const ublas::vector<T> &, const ublas::vector<T> &,
    ublas::vector<T> &, Method, const ublas::matrix<T> &, ublas::matrix<T> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void step_y(const T, const T, const T,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &, std::size_t, const T, std::size_t, const T,
    const T, std::size_t, const T, const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &,
    const ublas::vector<T> &, const ublas::vector<T> &, const ublas::vector<T> &,
    ublas::vector<T> &, Method, const ublas::matrix<T> &, ublas::matrix<T> &);

template<typename T, typename>
ublas::matrix<T> splitting_method(const T a, const T b, const T mu,
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
    const Method method, const std::function<void (const ublas::matrix<T> &)> &get_error)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h_1 = l_1 / n_1, h_2 = l_2 / n_2, tau = t / k_upper,
        sigma_a = a * tau / (h_1 * h_1), sigma_b = b * tau / (h_2 * h_2);
    if (n_1 < 2U || n_2 < 2U || !k_upper ||
        h_1 < EPSILON || h_1 < EPSILON || tau < EPSILON || !is_enum_includes(method))
    {
        throw std::logic_error("n_1 < 2 || n_2 < 2 || !k || "
            "h_1 < epsilon || h_2 < epsilon || tau < epsilon");
    }

    ublas::vector<T> d_i(n_1 + 1), d_j(n_2 + 1),
        a_i(n_1 + 1, sigma_a / (method == ALTERNATING_DIRECTION ? 2.0 : 1.0)),
        b_i(n_1 + 1, -(method == ALTERNATING_DIRECTION ? 1.0 : 2.0) * sigma_a - 1.0),
        c_i(n_1 + 1, sigma_a / (method == ALTERNATING_DIRECTION ? 2.0 : 1.0)),
        a_j(n_2 + 1, sigma_b / (method == ALTERNATING_DIRECTION ? 2.0 : 1.0)),
        b_j(n_2 + 1, -(method == ALTERNATING_DIRECTION ? 1.0 : 2.0) * sigma_b - 1.0),
        c_j(n_2 + 1, sigma_b / (method == ALTERNATING_DIRECTION ? 2.0 : 1.0));
    a_i(0) = c_i(n_1) = a_j(0) = c_j(n_2) = 0.0;
    b_i(0) = beta_3 - alpha_3 / h_1;
    c_i(0) = alpha_3 / h_1;
    a_i(n_1) = -alpha_4 / h_1;
    b_i(n_1) = beta_4 + alpha_4 / h_1;
    b_j(0) = beta_1 - alpha_1 / h_2;
    c_j(0) = alpha_1 / h_2;
    a_j(n_2) = -alpha_2 / h_2;
    b_j(n_2) = beta_2 + alpha_2 / h_2;

    ublas::matrix<T> u_k(n_1 + 1, n_2 + 1), u_k_minus_1_divides_2(n_1 + 1, n_2 + 1),
        u_k_minus_1(n_1 + 1, n_2 + 1);
    initial_approximation(a, b, mu, n_1, h_1, n_2, h_2, psi_x_y, u_k_minus_1);
    get_error(u_k_minus_1);
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        step_x(a, b, mu, f_x_y_t,
            n_1, h_1, n_2, h_2, tau, k, sigma_b, alpha_1, beta_1, alpha_2, beta_2,
            phi_1_x_t, phi_2_x_t, phi_3_y_t, phi_4_y_t, a_i, b_i, c_i, d_i,
            method, u_k_minus_1, u_k_minus_1_divides_2);
        step_y(a, b, mu, f_x_y_t,
            n_1, h_1, n_2, h_2, tau, k, sigma_a, alpha_3, beta_3, alpha_4, beta_4,
            phi_1_x_t, phi_2_x_t, phi_3_y_t, phi_4_y_t, a_j, b_j, c_j, d_j,
            method, u_k_minus_1_divides_2, u_k);
        get_error(u_k);
        u_k.swap(u_k_minus_1);
    }

    return u_k_minus_1;
}

static constexpr bool is_enum_includes(const Method method) noexcept
{
    switch (method)
    {
        case ALTERNATING_DIRECTION:
        case FRACTIONAL_STEP:
            return true;
        default:
            return false;
    }
}

template<typename T, typename>
static void initial_approximation(const T a, const T b, const T mu,
    const std::size_t n_1, const T h_1, const std::size_t n_2, const T h_2,
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
static void step_x(const T a, const T b, const T mu,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &f_x_y_t,
    const std::size_t n_1, const T h_1, const std::size_t n_2, const T h_2,
    const T tau, const std::size_t k, const T sigma_b,
    const T alpha_1, const T beta_1, const T alpha_2, const T beta_2,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_1_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_2_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_3_y_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_4_y_t,
    const ublas::vector<T> &a_i, const ublas::vector<T> &b_i, const ublas::vector<T> &c_i,
    ublas::vector<T> &d_i, const Method method, const ublas::matrix<T> &u_k_minus_1,
    ublas::matrix<T> &u_k_minus_1_divides_2)
{
    for (std::size_t j = 1; j < n_2; ++j)
    {
        d_i(0) = phi_3_y_t(a, b, mu, j * h_2, (k - 0.5) * tau);
        d_i(n_1) = phi_4_y_t(a, b, mu, j * h_2, (k - 0.5) * tau);
        for (std::size_t i = 1; i < n_1; ++i)
        {
            d_i(i) = -u_k_minus_1(i, j)
                - tau / 2.0 * f_x_y_t(a, b, mu, i * h_1, j * h_2, (k - 0.5) * tau);
            if (method == ALTERNATING_DIRECTION)
            {
                d_i(i) += sigma_b * u_k_minus_1(i, j)
                    - sigma_b / 2.0 * u_k_minus_1(i, j - 1)
                    - sigma_b / 2.0 * u_k_minus_1(i, j + 1);
            }
        }
        ublas::column(u_k_minus_1_divides_2, j) = thomas_algorithm(a_i, b_i, c_i, d_i);
    }
    for (std::size_t i = 0; i <= n_1; ++i)
    {
        u_k_minus_1_divides_2(i, 0) = (-alpha_1 * u_k_minus_1_divides_2(i, 1)
                + h_2 * phi_1_x_t(a, b, mu, i * h_1, (k - 0.5) * tau))
            / (beta_1 * h_2 - alpha_1);
        u_k_minus_1_divides_2(i, n_2) = (alpha_2 * u_k_minus_1_divides_2(i, n_2 - 1)
                + h_2 * phi_2_x_t(a, b, mu, i * h_1, (k - 0.5) * tau))
            / (beta_2 * h_2 + alpha_2);
    }
}

template<typename T, typename>
static void step_y(const T a, const T b, const T mu,
    const std::function<T (const T &, const T &, const T &,
        const T &, const T &, const T &)> &f_x_y_t,
    const std::size_t n_1, const T h_1, const std::size_t n_2, const T h_2,
    const T tau, const std::size_t k, const T sigma_a,
    const T alpha_3, const T beta_3, const T alpha_4, const T beta_4,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_1_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_2_x_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_3_y_t,
    const std::function<T (const T &, const T &, const T &, const T &, const T &)> &phi_4_y_t,
    const ublas::vector<T> &a_j, const ublas::vector<T> &b_j, const ublas::vector<T> &c_j,
    ublas::vector<T> &d_j, const Method method, const ublas::matrix<T> &u_k_minus_1_divides_2,
    ublas::matrix<T> &u_k)
{
    for (std::size_t i = 1; i < n_1; ++i)
    {
        d_j(0) = phi_1_x_t(a, b, mu, i * h_1, k * tau);
        d_j(n_2) = phi_2_x_t(a, b, mu, i * h_1, k * tau);
        for (std::size_t j = 1; j < n_2; ++j)
        {
            d_j(j) = -u_k_minus_1_divides_2(i, j)
                - tau / 2.0 * f_x_y_t(a, b, mu, i * h_1, j * h_2, k * tau);
            if (method == ALTERNATING_DIRECTION)
            {
                d_j(j) += sigma_a * u_k_minus_1_divides_2(i, j)
                    - sigma_a / 2.0 * u_k_minus_1_divides_2(i - 1, j)
                    - sigma_a / 2.0 * u_k_minus_1_divides_2(i + 1, j);
            }
        }
        ublas::row(u_k, i) = thomas_algorithm(a_j, b_j, c_j, d_j);
    }
    for (std::size_t j = 0; j <= n_2; ++j)
    {
        u_k(0, j) = (-alpha_3 * u_k(1, j)
            + h_1 * phi_3_y_t(a, b, mu, j * h_2, k * tau)) / (beta_3 * h_1 - alpha_3);
        u_k(n_1, j) = (alpha_4 * u_k(n_1 - 1, j)
            + h_1 * phi_4_y_t(a, b, mu, j * h_2, k * tau)) / (beta_4 * h_1 + alpha_4);
    }
}

#endif
