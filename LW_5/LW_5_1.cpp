#include <cmath>

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <boost/numeric/ublas/io.hpp>

//#include "parabolic_pde.hpp"
#include "../generic/lup.hpp"
#include "../generic/qr.hpp"
#include "../generic/thomas_algorithm.hpp"

typedef long double ldbl;

static constexpr ldbl PI_LDBL = 3.1415926535897932384626L;

constexpr ldbl u_0_t(const ldbl &) noexcept;
constexpr ldbl u_l_t(const ldbl &) noexcept;
constexpr ldbl u_x_0(const ldbl &) noexcept;
constexpr ldbl u_exact(const ldbl &, const ldbl &, const ldbl &) noexcept;

int main()
{
    try
    {
/*
        const ldbl a = 1.0, l = 1.0, t = 2.0;
        const std::size_t n = 4096U, k = 64U;
        const ublas::vector<ldbl>
            explicit_fdm_w = explicit_fdm<ldbl>(a, l, n, t, k, u_0_t, u_l_t, u_x_0),
            implicit_fdm_w = implicit_fdm<ldbl>(a, l, n, t, k, u_0_t, u_l_t, u_x_0);
        Удалить mingw из path
        for (std::size_t i = 0; i < explicit_fdm_w.size(); ++i)
        {
            std::cout << std::abs(explicit_fdm_w[i] - implicit_fdm_w[i]) << ' ';
        }
        std::cout << '\n';
*/
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

constexpr ldbl u_0_t(const ldbl &) noexcept
{
    return 0.0;
}

constexpr ldbl u_l_t(const ldbl &) noexcept
{
    return 0.0;
}

constexpr ldbl u_x_0(const ldbl &x) noexcept
{
    return std::sin(2.0 * PI_LDBL * x);
}

constexpr ldbl u_exact(const ldbl &a, const ldbl &x, const ldbl &t) noexcept
{
    return std::exp(-4.0 * PI_LDBL * PI_LDBL * a * t) * std::sin(2.0 * PI_LDBL * x);
}
