#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <cstddef>

#include <array>
#include <type_traits>
#include <vector>
#include <utility>

typedef unsigned uint;

static constexpr const uint COUNT = 3;

enum Method : uint
{
    RECTANGLE,
    TRAPEZOIDAL,
    SIMSON
};

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T rectangle_method(const std::vector<T> &, T (const T &), const T &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T trapezoidal_rule(const std::vector<T> &, const T &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T simpson_s_rule(const std::vector<T> &, const T &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
std::pair<T, T> runge_rule(const T &, const T &, const T &, const T &, uint) noexcept;

template<typename T, typename>
T rectangle_method(const std::vector<T> &x, T func(const T &), const T &step)
{
    T result = 0.0;
    for (std::size_t size = x.size(), i = 1; i < size; ++i)
    {
        result += func((x[i - 1] + x[i]) / 2.0) * step;
    }

    return result;
}

template<typename T, typename>
T trapezoidal_rule(const std::vector<T> &f_of_x, const T &step)
{
    T result = 0.0;
    for (std::size_t size = f_of_x.size(), i = 1; i < size; ++i)
    {
        result += f_of_x[i - 1] + f_of_x[i];
    }
    return result * step / 2.0;
}

template<typename T, typename>
T simpson_s_rule(const std::vector<T> &f_of_x, const T &step)
{
    static constexpr const size_t TWO = 2;
    const std::array<T, TWO> c = { 2., 4. };
    T result = f_of_x.front() + f_of_x.back();
    for (std::size_t size_minus_1 = f_of_x.size() - 1, i = 1; i < size_minus_1; ++i)
    {
        result += c[i % TWO] * f_of_x[i];
    }

    return result * step / 3.0;
}

template<typename T, typename>
std::pair<T, T> runge_rule(const T &f_hk, const T &hk, const T &f_h, const T &h,
    const uint method) noexcept
{

    static constexpr const std::array<T, COUNT> p = { 1.0, 2.0, 4.0 };
    const T k = hk / h,
        result = f_h + (f_h - f_hk) / (std::pow(k, p[method]) - 1.0);
    return std::pair<T, T>(result, std::pow(h, p[method] + 1.0));
}

#endif
