#include <cmath>

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <limits>

#include "../generic/header.hpp"
#include "parabolic_pde.hpp"

static constexpr ldbl u_0_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_l_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_x_0(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &, const ldbl &) noexcept;

static constexpr ldbl L = PI_LDBL;

int main()
{
    try
    {
        ldbl a, t;
        std::size_t n, k;
        std::cin >> a >> t >> n >> k;
        const ublas::vector<ldbl>
            explicit_fdm_w = explicit_fdm<ldbl>(a, L, n, t, k, u_0_t, u_l_t, u_x_0),
            implicit_fdm_w = implicit_fdm<ldbl>(a, L, n, t, k, u_0_t, u_l_t, u_x_0),
            crank_nicolson_w = crank_nicolson<ldbl>(a, L, n, t, k, u_0_t, u_l_t, u_x_0);
        const ldbl h = L / n;
        for (std::size_t i = 0; i < explicit_fdm_w.size(); ++i)
        {
            std::cout << std::fixed << std::setprecision(9) <<
                explicit_fdm_w[i] << ' ' << implicit_fdm_w[i] << ' ' <<
                crank_nicolson_w[i] << ' ' << u_exact(a, i * h, t) << '\n';
        }
        std::cout << '\n';
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

static constexpr ldbl u_0_t(const ldbl &a, const ldbl &t) noexcept
{
    return std::exp(-a * t);
}

static constexpr ldbl u_l_t(const ldbl &a, const ldbl &t) noexcept
{
    return -std::exp(-a * t);
}

static constexpr ldbl u_x_0(const ldbl &, const ldbl &x) noexcept
{
    return std::cos(x);
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &x, const ldbl &t) noexcept
{
    return std::exp(-a * t) * std::cos(x);
}
