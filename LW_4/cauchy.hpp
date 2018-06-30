#ifndef __CAUCHY_HPP__
#define __CAUCHY_HPP__

#include <cstddef>

#include <type_traits>
#include <vector>

typedef unsigned uint;

static constexpr std::size_t ADAMS_FACTOR_COUNT = 4;

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
void euler_method(T (const T &, const T &, const T &),
    const std::vector<T> &, const T &, const T &, const T &,
    std::vector<T> &, std::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
void mod_euler_method(T (const T &, const T &, const T &),
    const std::vector<T> &, const T &, const T &, const T &,
    std::vector<T> &, std::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
void runge_kutt_method(T (const T &, const T &, const T &),
    const std::vector<T> &, const T &, const T &, const T &,
    std::vector<T> &, std::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
void adams_method(T (const T &, const T &, const T &),
    const std::vector<T> &, const T &, const T &, const T &,
    std::vector<T> &, std::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T adams_itearation(const std::vector<T> &, std::size_t);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T adams_itearation(T (const T &, const T &, const T &),
    const std::vector<T> &, const std::vector<T> &, const std::vector<T> &,
    std::size_t);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
std::vector<T> runge_rule(const std::vector<T> &, const T &,
    const std::vector<T> &, const T &, uint);

template<typename T, typename>
void euler_method(T func(const T &, const T &, const T &),
    const std::vector<T> &x, const T &step, const T &y_begin, const T &dy_begin,
    std::vector<T> &y, std::vector<T> &dy)
{
    const std::size_t size = x.size();
    y.resize(size);
    dy.resize(size);
    y.front() = y_begin;
    dy.front() = dy_begin;
    for (std::size_t i = 1; i < size; ++i)
    {
        y[i] = y[i - 1] + step * dy[i - 1];
        dy[i] = dy[i - 1] + step * func(x[i - 1], y[i - 1], dy[i - 1]);
    }
}

template<typename T, typename>
void mod_euler_method(T func(const T &, const T &, const T &),
    const std::vector<T> &x, const T &step, const T &y_begin, const T &dy_begin,
    std::vector<T> &y, std::vector<T> &dy)
{
    const std::size_t size = x.size();
    y.resize(size);
    dy.resize(size);
    y.front() = y_begin;
    dy.front() = dy_begin;
    const T step_half = step / 2.;
    T x_half, y_half, dy_half;

    for (std::size_t i = 1; i < size; ++i)
    {
        x_half = x[i - 1] + step_half;
        y_half = y[i - 1] + step_half * dy[i - 1];
        dy_half = dy[i - 1] + step_half * func(x[i - 1], y[i - 1], dy[i - 1]);

        y[i] = y[i - 1] + step * dy_half;
        dy[i] = dy[i - 1] + step * func(x_half, y_half, dy_half);
    }
}

template<typename T, typename>
void runge_kutt_method(T func(const T &, const T &, const T &),
    const std::vector<T> &x, const T &step, const T &y_begin, const T &dy_begin,
    std::vector<T> &y, std::vector<T> &dy)
{
    const std::size_t size = x.size();
    y.resize(size);
    dy.resize(size);
    y.front() = y_begin;
    dy.front() = dy_begin;

    for (std::size_t i = 1; i < size; ++i)
    {
        const T
            k1 = step * dy[i - 1],
            l1 = step * func(x[i - 1], y[i - 1], dy[i - 1]),
            k2 = step * (dy[i - 1] + l1 / 2.),
            l2 = step * func(x[i - 1] + step / 2., y[i - 1] + k1 / 2., dy[i - 1] + l1 / 2.),
            k3 = step * (dy[i - 1] + l2 / 2.),
            l3 = step * func(x[i - 1] + step / 2., y[i - 1] + k2 / 2.,  dy[i - 1] + l2 / 2.),
            k4 = step * (dy[i - 1] + l3),
            l4 = step * func(x[i - 1] + step, y[i - 1] + k3, dy[i - 1] + l3);

        y[i] = y[i - 1] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
        dy[i] = dy[i - 1] + (l1 + 2. * l2 + 2. * l3 + l4) / 6.;
    }
}

template<typename T, typename>
void adams_method(T func(const T &, const T &, const T &),
    const std::vector<T> &x, const T &step, const T &y_begin, const T &dy_begin,
    std::vector<T> &y, std::vector<T> &dy)
{
    if (y.size() < ADAMS_FACTOR_COUNT || dy.size() < ADAMS_FACTOR_COUNT)
    {
        throw std::logic_error("y.size() < 4 || dy.size() < 4");
    }

    const std::size_t size = x.size();
    y.resize(size);
    dy.resize(size);
    y.front() = y_begin;
    dy.front() = dy_begin;

    for (std::size_t i = ADAMS_FACTOR_COUNT; i < size; ++i)
    {
        y[i] = y[i - 1] + step * adams_itearation(dy, i);
        dy[i] = dy[i - 1] + step * adams_itearation(func, x, y, dy, i);
    }
}

template<typename T, typename>
std::vector<T> runge_rule(const std::vector<T> &f_hk, const T &hk,
    const std::vector<T> &f_h, const T &h, const uint p)
{
    const std::size_t size = f_hk.size();
    std::vector<T> result(size);
    const T k = std::round(hk / h);
    for (std::size_t i = 0; i < size; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(k * i);
        result[i] = f_h[idx] + (f_h[idx] - f_hk[i]) / (std::pow(k, p) - 1.);
    }
    return result;
}

template<typename T, typename>
T adams_itearation(const std::vector<T> &value, const std::size_t i)
{
    return (55. * value[i - 1] -
            59. * value[i - 2] +
            37. * value[i - 3] -
             9. * value[i - 4]) / 24.;
}

template<typename T, typename>
T adams_itearation(T func(const T &, const T &, const T &),
    const std::vector<T> &x, const std::vector<T> &y, const std::vector<T> &dy,
    const std::size_t i)
{
    return (55. * func(x[i - 1], y[i - 1], dy[i - 1]) -
            59. * func(x[i - 2], y[i - 2], dy[i - 2]) +
            37. * func(x[i - 3], y[i - 3], dy[i - 3]) -
             9. * func(x[i - 4], y[i - 4], dy[i - 4])) / 24.;
}

#endif
