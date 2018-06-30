#include <cmath>

#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "cauchy.hpp"

typedef double dbl;

constexpr dbl func_f(const dbl &, const dbl &, const dbl &) noexcept;
constexpr dbl exact_func_f(const dbl &) noexcept;

template <typename T>
std::ostream &out_vector(std::ostream &, std::size_t, std::size_t,
    const std::vector<T> &);

template <typename T>
void cauchy(std::ostream &, T (const T &, const T &, const T &),
    const T &, const T &, const T &, const T &, const T &, T (const T &));

int main()
{
    dbl begin, end, h, y_begin, dy_begin;
    std::cin >> begin >> end >> h >> y_begin >> dy_begin;
    try
    {
        cauchy(std::cout, func_f, begin, end, h, y_begin, dy_begin, exact_func_f);
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }
    std::cout << std::endl;

    return 0;
}

template <typename T>
void cauchy(std::ostream &os, T func(const T &, const T &, const T &),
    const T &begin, const T &end, const T &step,
    const T &y_begin, const T &dy_begin, T exact_func(const T &))
{
    static constexpr std::size_t WIDTH = 1, PRECISION = 6;
    const std::size_t size = static_cast<std::size_t>((end - begin) / step + 1.);
    if (begin >= end || step < std::numeric_limits<T>::epsilon() ||
        step > end - begin || size < ADAMS_FACTOR_COUNT)
    {
        throw std::logic_error("begin >= end || step < epsilon || point_count < 4");
    }

    std::vector<T> x(size), exact_y(size);
    os << "[" << begin << ", " << end << "]" << '\n';
    os << "With step = " << step << '\n';

    for (std::size_t i = 0; i < size; ++i)
    {
        const T value = i * step;
        x[i] = value;
        exact_y[i] = exact_func(value);
    }
    os << "x" << '\n';
    out_vector(os, PRECISION, WIDTH, x);
    os << "exact_y" << '\n';
    out_vector(os, PRECISION, WIDTH, exact_y);

    std::vector<T> y, dy;

    euler_method(func, x, step, y_begin, dy_begin, y, dy);
    os << "euler_method" << '\n';
    out_vector(os, PRECISION, WIDTH, y);

    mod_euler_method(func, x, step, y_begin, dy_begin, y, dy);
    os << "mod_euler_method" << '\n';
    out_vector(os, PRECISION, WIDTH, y);

    runge_kutt_method(func, x, step, y_begin, dy_begin, y, dy);
    os << "runge_kutt_method" << '\n';
    out_vector(os, PRECISION, WIDTH, y);

    adams_method(func, x, step, y_begin, dy_begin, y, dy);
    os << "adams_method" << '\n';
    out_vector(os, PRECISION, WIDTH, y);
}

constexpr dbl func_f(const dbl &x, const dbl &y, const dbl &dy) noexcept
{
    return -4. * x * dy - (4. * x * x + 2.) * y;
}

constexpr dbl exact_func_f(const dbl &x) noexcept
{
    return (1. + x) * exp(-x * x);
}

template <typename T>
std::ostream &out_vector(std::ostream &os, const std::size_t precision, const std::size_t width,
    const std::vector<T> &value)
{
    const std::size_t old_precision = os.precision(), old_width = os.width();
    os << std::fixed << std::setprecision(precision) << std::setw(width);
    for (const T &v_i : value)
    {
        os << v_i << ' ';
    }
    os << std::setprecision(old_precision) << std::setw(old_width);
    return os << '\n' << std::endl;
}
