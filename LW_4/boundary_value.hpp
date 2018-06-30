#ifndef __BOUNDARY_VALUE_HPP__
#define __BOUNDARY_VALUE_HPP__

#include <cstddef>

#include <functional>
#include <type_traits>
#include <vector>

#include "cauchy.hpp"
#include "../generic/thomas_algorithm.hpp"

typedef unsigned uint;

static constexpr std::size_t THOMAS_ALGORITHM_COUNT = 3;

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
void shooting_method_first_type_cond(T (const T &, const T &, const T &),
    const std::vector<T> &, const T &, const T &, const T &,
    const T &, const T &, const T &, std::vector<T> &, std::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
void finite_diff_method_first_type_cond(T (const T &), T (const T &), T (const T &),
    const std::vector<T> &, const T &, const T &, const T &,
    std::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T secant_method(const std::function<T (const T &)> &,
    const T &, const T &, const T &);

template<typename T, typename>
void shooting_method_first_type_cond(T func(const T &, const T &, const T &),
    const std::vector<T> &x, const T &step, const T &y_begin, const T &y_end,
    const T &dy_begin_approx_0, const T &dy_begin_approx_1,
    const T &epsilon, std::vector<T> &y, std::vector<T> &dy)
{
    const std::function<T (const T &)> func_cauchy =
        [&](const T &dy_begin_value) -> T
        {
            runge_kutt_method(func, x, step, y_begin, dy_begin_value, y, dy);
            return y.back() - y_end;
        };
    secant_method(func_cauchy, dy_begin_approx_0, dy_begin_approx_1, epsilon);
}

template<typename T, typename>
void finite_diff_method_first_type_cond(T p(const T &), T q(const T &), T f(const T &),
    const std::vector<T> &x, const T &step, const T &y_begin, const T &y_end,
    std::vector<T> &y)
{
    const std::size_t size = x.size(), size_minus_2 = x.size() - 2;
    if (size < THOMAS_ALGORITHM_COUNT + 2)
    {
        throw std::logic_error("point_count < 5");
    }

    y.resize(size);
    ublas::vector<T> a(size_minus_2), b(size_minus_2), c(size_minus_2), d(size_minus_2);
    for (std::size_t i = 0; i < size_minus_2; ++i)
    {
        a[i] = 1. / (step * step) - p(x[i + 1]) / (2. * step);
        b[i] = -2. / (step * step) + q(x[i + 1]);
        c[i] = 1. / (step * step) + p(x[i + 1]) / (2. * step);
        d[i] = f(x[i + 1]);
    }

    d[0] -= a[0] * y_begin;
    d[size_minus_2 - 1] -= c[size_minus_2 - 1] * y_end;

    a[0] = 0.;
    c[size_minus_2 - 1] = 0.;

    ublas::vector<T> temp_y = thomas_algorithm(a, b, c, d);
    y.front() = y_begin;
    y.back() = y_end;
    std::copy(temp_y.cbegin(), temp_y.cend(), y.begin() + 1);
}

template<typename T, typename>
T secant_method(const std::function<T (const T &)> &func,
    const T &x_0, const T &x_1, const T &epsilon)
{
    T x_k_minus_1, x_k = x_0, x_k_plus_1 = x_1, f_x_k_minus_1, f_x_k = func(x_0);
    do
    {
        x_k_minus_1 = x_k;
        f_x_k_minus_1 = f_x_k;

        x_k = x_k_plus_1;
        f_x_k = func(x_k);

        x_k_plus_1 = x_k - f_x_k * (x_k - x_k_minus_1) / (f_x_k - f_x_k_minus_1);
    } while (std::fabs(x_k_plus_1 - x_k) >= epsilon);

    return x_k_plus_1;
}

#endif
