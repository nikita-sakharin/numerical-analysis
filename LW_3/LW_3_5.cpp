#include <cmath>

#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "integration.hpp"

typedef double dbl;

constexpr dbl func(const dbl &) noexcept;

template <typename T>
std::array<T, COUNT> integration(std::ostream &, T (const T &),
    const T &, const T &, const T &);

int main()
{
    dbl begin, end, h1, h2, exact_value;
    std::cin >> begin >> end >> h1 >> h2 >> exact_value;
    try
    {
        std::array<dbl, COUNT> f1, f2;
        f1 = integration(std::cout, func, begin, end, h1);
        f2 = integration(std::cout, func, begin, end, h2);

        std::pair<dbl, dbl> runge;
        runge = runge_rule(f1[RECTANGLE], h1, f2[RECTANGLE], h2, RECTANGLE);
        std::cout << "Runge correction for rectangle method = " << runge.first << '\n';
        std::cout << "Calculation error <= " <<
            std::fabs(runge.first - exact_value) << '\n';

        runge = runge_rule(f1[TRAPEZOIDAL], h1, f2[TRAPEZOIDAL], h2, TRAPEZOIDAL);
        std::cout << "Runge correction for trapezoidal rule = " << runge.first << '\n';
        std::cout << "Calculation error <= " <<
            std::fabs(runge.first - exact_value) << '\n';

        runge = runge_rule(f1[SIMSON], h1, f2[SIMSON], h2, SIMSON);
        std::cout << "Runge correction for simpson's rule   = " << runge.first << '\n';
        std::cout << "Calculation error <= " <<
            std::fabs(runge.first - exact_value) << '\n';

    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }
    std::cout << std::endl;

    return 0;
}

template <typename T>
std::array<T, COUNT> integration(std::ostream &os, T func(const T &),
    const T &begin, const T &end, const T &step)
{
    const std::size_t size = static_cast<std::size_t>((end - begin) / step + 1.);
    if (begin >= end || step < std::numeric_limits<T>::epsilon() || size < 2)
    {
        throw std::logic_error("begin >= end || step < epsilon || point_count < 2");
    }

    std::vector<T> x(size), f_of_x(size);
    os << "[" << begin << ", " << end << "]" << '\n';
    os << "With step step = " << step << '\n';

    T value = begin;
    for (std::size_t i = 0; i < size; ++i, value += step)
    {
        x[i] = value;
        f_of_x[i] = func(value);
    }

    std::array<T, COUNT> result = {
        rectangle_method(x, func, step),
        trapezoidal_rule(f_of_x, step),
        simpson_s_rule(f_of_x, step)
    };

    os << "Rectangle method: F = " << result[RECTANGLE] << '\n';
    os << "Trapezoidal rule: F = " << result[TRAPEZOIDAL] << '\n';
    os << "Simpson's rule:   F = " << result[SIMSON] << '\n';
    os << std::endl;
    return result;
}

constexpr dbl func(const dbl &x) noexcept
{
    return x /
        ((2.0 * x + 7.0) * (3.0 * x + 4.0));
}
