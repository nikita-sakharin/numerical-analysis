#ifndef __CUBIC_SPLINE_HPP__
#define __CUBIC_SPLINE_HPP__

#include <cstddef>

#include <algorithm>
#include <array>
#include <type_traits>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../generic/thomas_algorithm.hpp"

using namespace boost::numeric;

constexpr static const std::size_t ORDER = 3;
constexpr static const std::size_t ORDER_PLUS_1 = 4;

enum Coef : std::size_t
{
    A,
    B,
    C,
    D
};

template <typename T, std::size_t ORDER>
using Polynom = std::array<T, ORDER + 1>;

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T cubic_spline(std::ostream &, const ublas::vector<T> &, const ublas::vector<T> &, const T);

template<typename T>
ublas::vector<Polynom<T, ORDER>>
    coefficient(const ublas::vector<T> &, const ublas::vector<T> &);

template<typename T>
ublas::vector<T> compute_h(const ublas::vector<T> &);

template<typename T>
void three_diagonal_matrix(const ublas::vector<T> &, const ublas::vector<T> &,
    ublas::vector<T> &, ublas::vector<T> &, ublas::vector<T> &, ublas::vector<T> &);

template<typename T>
T interpolate(const ublas::vector<T> &, const T,
    const ublas::vector<Polynom<T, ORDER>> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
std::ostream &polynom_interal_out(std::ostream &, const ublas::vector<T> &,
    const ublas::vector<Polynom<T, ORDER>> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T cubic_spline(std::ostream &os, const ublas::vector<T> &x, const ublas::vector<T> &f_of_x,
    const T x_star)
{
    const std::size_t size = x.size();
    if (size != f_of_x.size() || size < 5 || x_star < x[0])
    {
        throw std::domain_error("Can't be solved");
    }

    const ublas::vector<Polynom<T, ORDER>> p =
        std::move(coefficient(x, f_of_x));
    polynom_interal_out(os, x, p);
    return interpolate(x, x_star, p);
}

template<typename T>
ublas::vector<Polynom<T, ORDER>>
    coefficient(const ublas::vector<T> &x, const ublas::vector<T> &f_of_x)
{
    const ublas::vector<T> h = std::move(compute_h(x));
    ublas::vector<T> a, b, c, d;
    three_diagonal_matrix(f_of_x, h, a, b, c, d);
    ublas::vector<T> p_c = std::move(thomas_algorithm(a, b, c, d));
    const std::size_t size = p_c.size() + 1;
    ublas::vector<Polynom<T, ORDER>> p(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        p[i][C] = i ? p_c[i - 1] : 0.;
    }
    for (std::size_t i = 0; i < size; ++i)
    {
        p[i][A] = f_of_x[i];
        if (i + 1 < size)
        {
            p[i][B] = (f_of_x[i + 1] - f_of_x[i]) / h[i + 1] -
                h[i + 1] * (p[i + 1][C] + 2. * p[i][C]) / 3.;
            p[i][D] = (p[i + 1][C] - p[i][C])
                / (3.0 * h[i + 1]);
        }
        else
        {
            p[i][B] = (f_of_x[i + 1] - f_of_x[i]) / h[i + 1] -
                2. / 3. * h[i + 1] * p[i][C];
            p[i][D] = -p[i][C] / (3.0 * h[i + 1]);
        }
    }

    return p;
}

template<typename T>
ublas::vector<T> compute_h(const ublas::vector<T> &x)
{
    const std::size_t size = x.size();
    ublas::vector<T> h(size);
    h[0] = std::numeric_limits<T>::quiet_NaN();
    for (std::size_t i = 1; i < size; ++i)
    {
        if (x[i] <= x[i - 1])
        {
            throw std::domain_error("Unsorted elements or two equal");
        }
        h[i] = x[i] - x[i - 1];
    }

    return h;
}

template<typename T>
void three_diagonal_matrix(const ublas::vector<T> &f_of_x, const ublas::vector<T> &h,
    ublas::vector<T> &a, ublas::vector<T> &b, ublas::vector<T> &c,
    ublas::vector<T> &d)
{
    const std::size_t size = f_of_x.size() - 2;
    a.resize(size);
    b.resize(size);
    c.resize(size);
    d.resize(size);

    for (std::size_t i = 0; i < size; ++i)
    {
        if (i)
        {
            a[i] = h[i + 1];
        }
        b[i] = 2. * (h[i + 1] + h[i + 2]);
        if (i + 1 < size)
        {
            c[i] = h[i + 2];
        }
        d[i] = 3. *
            ((f_of_x[i + 2] - f_of_x[i + 1]) / h[i + 2] -
            (f_of_x[i + 1] - f_of_x[i]) / h[i + 1]);
    }
}

template<typename T>
T interpolate(const ublas::vector<T> &x, const T x_star,
    const ublas::vector<Polynom<T, ORDER>> &p)
{
    std::size_t i = std::lower_bound(x.cbegin(), x.cend(), x_star) - x.cbegin();
    if (i == x.size() || x_star != x[i])
    {
        --i;
    }
    const T delta = x_star - x[i];
    T prod = 1., result = 0.;
    for (std::size_t j = 0; j < ORDER_PLUS_1; ++j)
    {
        result += p[i][j] * prod;
        prod *= delta;
    }

    return result;
}


template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
std::ostream &polynom_interal_out(std::ostream &os, const ublas::vector<T> &x,
    const ublas::vector<Polynom<T, ORDER>> &p)
{
    const std::size_t size = x.size();
    for (std::size_t i = 0; i < size - 1; ++i)
    {

        os << "[" << x(i) << ", " << x(i + 1) << "] : "
            << p[i][A] << " + "
            << p[i][B] << " * x + "
            << p[i][C] << " * x^2 + "
            << p[i][D] << " * x^3\n";
    }
    return os << '\n';
}

#endif
