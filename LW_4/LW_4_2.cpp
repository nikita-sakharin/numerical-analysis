#include <cmath>

#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "boundary_value.hpp"

typedef long double ldbl;

constexpr ldbl func_F(const ldbl &, const ldbl &, const ldbl &) noexcept;

constexpr ldbl func_F_p(const ldbl &) noexcept;
constexpr ldbl func_F_q(const ldbl &) noexcept;
constexpr ldbl func_F_f(const ldbl &) noexcept;

constexpr ldbl exact_func_F(const ldbl &) noexcept;

template <typename T>
std::ostream &out_vector(std::ostream &, std::size_t, std::size_t,
    const std::vector<T> &);

template <typename T>
void boundary_value_first_type_cond(std::ostream &,
    T (const T &, const T &, const T &),
    T (const T &), T (const T &), T (const T &), const T &, const T &,
    const T &, const T &, const T &, T (const T &));

int main()
{
    ldbl begin, end, h, y_begin, y_end;
    std::cin >> begin >> end >> h >> y_begin >> y_end;
    try
    {
        boundary_value_first_type_cond(std::cout, func_F,
            func_F_p, func_F_q, func_F_f, begin, end, h,
            y_begin, y_end, exact_func_F);
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }
    std::cout << std::endl;
    return 0;
}

template <typename T>
void boundary_value_first_type_cond(std::ostream &os,
    T func(const T &, const T &, const T &),
    T func_p(const T &), T func_q(const T &), T func_f(const T &),
    const T &begin, const T &end,
    const T &step, const T &y_begin, const T &y_end, T exact_func(const T &))
{
    static constexpr std::size_t WIDTH = 1, PRECISION = 6;
    static constexpr T EPSILON = 1E-9;

    const std::size_t size = static_cast<std::size_t>((end - begin) / step + 1.);
    if (begin >= end || step < std::numeric_limits<T>::epsilon() || step > end - begin || size < 5)
    {
        throw std::logic_error("begin >= end || step < epsilon || point_count < 5");
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
    const T dy_begin_approx_0 = (y_end - y_begin) / (end - begin),
        dy_begin_approx_1 = (dy_begin_approx_0 + exact_func(end)) / 2.;

    shooting_method_first_type_cond(func, x, step, y_begin, y_end,
        dy_begin_approx_0, dy_begin_approx_1, EPSILON, y, dy);
    os << "shooting_method" << '\n';
    out_vector(os, PRECISION, WIDTH, y);

    finite_diff_method_first_type_cond(func_p, func_q, func_f, x, step, y_begin, y_end, y);
    os << "finite_diff_method_first_type_cond" << '\n';
    out_vector(os, PRECISION, WIDTH, y);
}

constexpr ldbl func_F(const ldbl &x, const ldbl &y, const ldbl &) noexcept
{
    return 2. * (1. + std::tan(x) * std::tan(x)) * y;
}

constexpr ldbl exact_func_F(const ldbl &x) noexcept
{
    return -std::tan(x);
}

constexpr ldbl func_F_p(const ldbl &) noexcept
{
    return 0.;
}

constexpr ldbl func_F_q(const ldbl &x) noexcept
{
    return -2. * (1. + std::tan(x) * std::tan(x));
}

constexpr ldbl func_F_f(const ldbl &) noexcept
{
    return 0.;
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
