#include <cmath>

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../generic/header.hpp"
#include "parabolic_pde.hpp"

static constexpr ldbl u_0_t(const ldbl &, const ldbl &, const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_l_t(const ldbl &, const ldbl &, const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_x_0(const ldbl &, const ldbl &, const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &, const ldbl &) noexcept;

static constexpr ldbl L = PI_LDBL;
static constexpr ldbl ALPHA = 1.0,
                      BETA  = 1.0,
                      GAMMA = 1.0,
                      DELTA = 1.0;
/*
static constexpr ldbl L = PI_LDBL;
static constexpr ldbl ALPHA = 0.0,
                      BETA  = 1.0,
                      GAMMA = 0.0,
                      DELTA = 1.0;
*/
int main()
{
    try
    {
        ldbl a, b, c, t;
        std::size_t n, k;
        std::cin >> a >> b >> c >> t >> n >> k;
        const ublas::vector<ldbl>
            explicit_fdm_w = explicit_fdm<ldbl>(a, b, c,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA, u_0_t, u_l_t, u_x_0),
            implicit_fdm_w = implicit_fdm<ldbl>(a, b, c,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA, u_0_t, u_l_t, u_x_0),
            crank_nicolson_w = crank_nicolson<ldbl>(a, b, c,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA, u_0_t, u_l_t, u_x_0);
        const ldbl h = L / n;
        for (std::size_t i = 0; i < explicit_fdm_w.size(); ++i)
        {
            std::cout << std::fixed << std::setprecision(18) << i * h << ',' <<
                u_exact(a, b, c, i * h, t) << ',' << explicit_fdm_w[i] << ',' <<
                implicit_fdm_w[i] << ',' << crank_nicolson_w[i] << '\n';
        }
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

static constexpr ldbl u_0_t(const ldbl &a, const ldbl &b, const ldbl &c, const ldbl &t) noexcept
{
    return std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
}

static constexpr ldbl u_l_t(const ldbl &a, const ldbl &b, const ldbl &c, const ldbl &t) noexcept
{
    return -std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
}

static constexpr ldbl u_x_0(const ldbl &, const ldbl &, const ldbl &, const ldbl &x) noexcept
{
    return std::sin(x);
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &x, const ldbl &t) noexcept
{
    return std::exp((c - a) * t) * std::sin(x + b * t);
}
/*
static constexpr ldbl u_0_t(const ldbl &a, const ldbl &, const ldbl &, const ldbl &t) noexcept
{
    return std::exp(-a * t);
}

static constexpr ldbl u_l_t(const ldbl &a, const ldbl &, const ldbl &, const ldbl &t) noexcept
{
    return -std::exp(-a * t);
}

static constexpr ldbl u_x_0(const ldbl &, const ldbl &, const ldbl &, const ldbl &x) noexcept
{
    return std::cos(x);
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &, const ldbl &,
    const ldbl &x, const ldbl &t) noexcept
{
    return std::exp(-a * t) * std::cos(x);
}
*/
