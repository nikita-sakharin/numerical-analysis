#include <cmath>

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../generic/header.hpp"
#include "elliptic_pde.hpp"

static constexpr ldbl f_x_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl phi_1_y(const ldbl &) noexcept;
static constexpr ldbl phi_2_y(const ldbl &) noexcept;
static constexpr ldbl phi_3_x(const ldbl &) noexcept;
static constexpr ldbl phi_4_x(const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &) noexcept;

static constexpr ldbl L_1 = PI_2_LDBL,
                      L_2 = PI_2_LDBL;
static constexpr ldbl ALPHA_1 = 0.0,
                      BETA_1  = 1.0,
                      ALPHA_2 = 0.0,
                      BETA_2  = 1.0,
                      ALPHA_3 = 0.0,
                      BETA_3  = 1.0,
                      ALPHA_4 = 0.0,
                      BETA_4  = 1.0,

int main()
{
    try
    {
        ldbl a, b, c;
        std::size_t n,_1 n_2;
        std::cin >> a >> b >> c >> n_1 >> n_2;
        const ublas::vector<ldbl>
            liebmann_fdm_u = liebmann_fdm<ldbl>(a, b, c, f_x_t,
                L_1, n_1, L_2, n_2, ALPHA_1, BETA_1, ALPHA_2, BETA_2,
                ALPHA_3, BETA_3, ALPHA_4, BETA_4,
                phi_1_y, phi_2_y, phi_3_x, phi_4_x),
            seidel_fdm_u = seidel_fdm<ldbl>(a, b, c, f_x_t,
                L_1, n_1, L_2, n_2, ALPHA_1, BETA_1, ALPHA_2, BETA_2,
                ALPHA_3, BETA_3, ALPHA_4, BETA_4,
                phi_1_y, phi_2_y, phi_3_x, phi_4_x);
        const ldbl h = L / n;
        for (std::size_t i = 0; i < explicit_fdm_u.size(); ++i)
        {
            std::cout << std::fixed << std::setprecision(18) << i * h << ',' <<
                u_exact(a, b, c, i * h, t) << ',' << explicit_fdm_u[i] << ',' <<
                implicit_fdm_u[i] << '\n';
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

static constexpr ldbl phi_1_y(const ldbl &t) noexcept
{
    return ;
}

static constexpr ldbl phi_2_y(const ldbl &y) noexcept
{
    return ;
}

static constexpr ldbl phi_3_x(const ldbl &x) noexcept
{
    return ;
}

static constexpr ldbl phi_4_x(const ldbl &x) noexcept
{
    return ;
}

static constexpr ldbl u_exact(const ldbl &x, const ldbl &y) noexcept
{
    return std::exp(-x - y) * std::cos(x) * std::cos(y);
}
