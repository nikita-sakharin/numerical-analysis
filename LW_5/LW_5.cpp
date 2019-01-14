#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../generic/header.hpp"
#include "parabolic_pde.hpp"

static constexpr ldbl f_x_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl phi_0_t(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &) noexcept;
static constexpr ldbl phi_l_t(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &) noexcept;
static constexpr ldbl psi_x(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &, const ldbl &) noexcept;

static constexpr ldbl L = PI_LDBL;
static constexpr ldbl ALPHA = 1.0,
                      BETA  = 1.0,
                      GAMMA = 1.0,
                      DELTA = 1.0;

int main(int argc, const char *argv[])
{
    try
    {
        ldbl a, b, c, t;
        std::size_t n, k;
        uint boundary;
        std::cin >> a >> b >> c >> t >> n >> k >> boundary;

        const ldbl h = L / n, tau = t / k;
        std::vector<ldbl> explicit_fdm_error,
            implicit_fdm_error,
            crank_nicolson_error;
        const ublas::vector<ldbl>
            explicit_fdm_u = explicit_fdm<ldbl>(a, b, c, f_x_t,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA,
                phi_0_t, phi_l_t, psi_x, static_cast<NumDiff>(boundary),
                [&] (const ublas::vector<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size = u_k.size(),
                        j = explicit_fdm_error.size();
                    for (std::size_t i = 0; i < size; ++i)
                    {
                        error = std::max(error, std::abs(u_k[i] - u_exact(a, b, c,
                            i * h, j * tau)));
                    }
                    explicit_fdm_error.push_back(error);
                }),
            implicit_fdm_u = implicit_fdm<ldbl>(a, b, c, f_x_t,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA,
                phi_0_t, phi_l_t, psi_x, static_cast<NumDiff>(boundary),
                [&] (const ublas::vector<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size = u_k.size(),
                        j = implicit_fdm_error.size();
                    for (std::size_t i = 0; i < size; ++i)
                    {
                        error = std::max(error, std::abs(u_k[i] - u_exact(a, b, c,
                            i * h, j * tau)));
                    }
                    implicit_fdm_error.push_back(error);
                }),
            crank_nicolson_u = crank_nicolson<ldbl>(a, b, c, f_x_t,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA,
                phi_0_t, phi_l_t, psi_x, static_cast<NumDiff>(boundary),
                [&] (const ublas::vector<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size = u_k.size(),
                        j = crank_nicolson_error.size();
                    for (std::size_t i = 0; i < size; ++i)
                    {
                        error = std::max(error, std::abs(u_k[i] - u_exact(a, b, c,
                            i * h, j * tau)));
                    }
                    crank_nicolson_error.push_back(error);
                });

        for (std::size_t i = 0; i < explicit_fdm_u.size(); ++i)
        {
            std::cout << std::fixed << std::setprecision(18) << i * h << ',' <<
                u_exact(a, b, c, i * h, t) << ',' << explicit_fdm_u[i] << ',' <<
                implicit_fdm_u[i] << ',' << crank_nicolson_u[i] << '\n';
        }
        if (argc > 1)
        {
            std::fstream error_stream(argv[1], std::ios_base::out | std::ios_base::trunc);
            for (std::size_t i = 0; i < explicit_fdm_error.size(); ++i)
            {
                error_stream << std::fixed << std::setprecision(18) << i * tau << ',' <<
                    explicit_fdm_error[i] << ',' << implicit_fdm_error[i] << ',' <<
                    crank_nicolson_error[i] << '\n';
            }
        }
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

static constexpr ldbl f_x_t(const ldbl &, const ldbl &) noexcept
{
    return 0.0;
}

static constexpr ldbl phi_0_t(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &t) noexcept
{
    return std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
}

static constexpr ldbl phi_l_t(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &t) noexcept
{
    return -std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
}

static constexpr ldbl psi_x(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &x) noexcept
{
    return std::sin(x);
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &x, const ldbl &t) noexcept
{
    return std::exp((c - a) * t) * std::sin(x + b * t);
}
