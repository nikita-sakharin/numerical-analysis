#ifndef __PARABOLIC_PDE_HPP__
#define __PARABOLIC_PDE_HPP__

#include <cstddef>

#include <type_traits>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../generic/thomas_algorithm.hpp"

template <typename T>
using Vector = boost::numeric::ublas::vector<T>;

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
Vector<T> explicit_fdm(const T, const T, std::size_t,
    const T, std::size_t, const std::function<T (const T &)> &,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
Vector<T> implicit_fdm(const T, const T, std::size_t,
    const T, std::size_t, const std::function<T (const T &)> &,
    const std::function<T (const T &)> &, const std::function<T (const T &)> &);

template<typename T, typename>
Vector<T> explicit_fdm(const T a, const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const std::function<T (const T &)> &u_0_t,
    const std::function<T (const T &)> &u_l_t,
    const std::function<T (const T &)> &u_x_0)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, a_sq = a * a;
    if (!n_upper || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("!n || !k || h < epsilon || tau < epsilon");
    }

    Vector<Vector<T>> w_h_tau(k_upper + 1, Vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = u_x_0(j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        w_h_tau[k][0] = u_0_t(k * tau);
        w_h_tau[k][n_upper] = u_l_t(k * tau);
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            w_h_tau[k][j] = w_h_tau[k - 1][j] + a_sq * tau *
                (w_h_tau[k - 1][j + 1] -
                2.0 * w_h_tau[k - 1][j] +
                w_h_tau[k - 1][j - 1]);
        }
    }

    return w_h_tau[k_upper];
}

template<typename T, typename>
Vector<T> implicit_fdm(const T a, const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const std::function<T (const T &)> &u_0_t,
    const std::function<T (const T &)> &u_l_t,
    const std::function<T (const T &)> &u_x_0)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper,
        a_sq = a * a, sigma = a_sq * tau / (h * h);
    const std::size_t n_upper_minus_1 = n_upper - 1;
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    Vector<Vector<T>> w_h_tau(k_upper + 1, Vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = u_x_0(j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        Vector<T> d_j(n_upper_minus_1),
            a_j(n_upper_minus_1, sigma),
            b_j(n_upper_minus_1, -(1.0 + 2.0 * sigma)),
            c_j(n_upper_minus_1, sigma);
        a_j[0] = c_j[n_upper_minus_1 - 1] = 0.0;
        d_j[0] =
            -(w_h_tau[k - 1][1] + sigma * u_0_t(k * tau));
        d_j[n_upper_minus_1 - 1] =
            -(w_h_tau[k - 1][n_upper - 1] + sigma * u_l_t(k * tau));
        for (std::size_t j = 1; j < n_upper_minus_1 - 1; ++j)
        {
            d_j[j] = -w_h_tau[k - 1][j + 1];
        }
        w_h_tau[k] = thomas_algorithm(a_j, b_j, c_j, d_j);
        w_h_tau[k].resize(n_upper + 1);
        std::move_backward(w_h_tau.begin(), w_h_tau.end() - 2, w_h_tau.end() - 1);
        w_h_tau[k][0] = u_0_t(k * tau);
        w_h_tau[k][n_upper] = u_l_t(k * tau);
    }

    return w_h_tau[k_upper];
}

template<typename T, typename>
Vector<T> crank_nicolson(const T a, const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const std::function<T (const T &)> &u_0_t,
    const std::function<T (const T &)> &u_l_t,
    const std::function<T (const T &)> &u_x_0)
{
    static constexpr T THETA = 0.5, EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper,
        a_sq = a * a, sigma = a_sq * tau / (h * h);
    const std::size_t n_upper_minus_1 = n_upper - 1;
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON)
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    Vector<Vector<T>> w_h_tau(k_upper + 1, Vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = u_x_0(j * h);
    }
    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        Vector<T> d_j(n_upper_minus_1),
            a_j(n_upper_minus_1, sigma),
            b_j(n_upper_minus_1, -(1.0 + 2.0 * sigma)),
            c_j(n_upper_minus_1, sigma);
        a_j[0] = c_j[n_upper_minus_1 - 1] = 0.0;
        d_j[0] =
            -(w_h_tau[k - 1][1] + sigma * u_0_t(k * tau));
        d_j[n_upper_minus_1 - 1] =
            -(w_h_tau[k - 1][n_upper - 1] + sigma * u_l_t(k * tau));
        for (std::size_t j = 1; j < n_upper_minus_1 - 1; ++j)
        {
            d_j[j] = -w_h_tau[k - 1][j + 1];
        }
        w_h_tau[k] = thomas_algorithm(a_j, b_j, c_j, d_j);
        w_h_tau[k].resize(n_upper + 1);
        std::move_backward(w_h_tau.begin(), w_h_tau.end() - 2, w_h_tau.end() - 1);
        w_h_tau[k][0] = u_0_t(k * tau);
        w_h_tau[k][n_upper] = u_l_t(k * tau);
    }

    return w_h_tau[k_upper];
}

#endif
