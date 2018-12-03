#include <cmath>

#include <fstream>
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
/*
static constexpr ldbl A = -2.0,
                      B = -2.0,
                      C = -4.0;
static constexpr ldbl L_1 = PI_2_LDBL,
                      L_2 = PI_2_LDBL;
static constexpr ldbl ALPHA_1 = 0.0,
                      BETA_1  = 1.0,
                      ALPHA_2 = 0.0,
                      BETA_2  = 1.0,
                      ALPHA_3 = 0.0,
                      BETA_3  = 1.0,
                      ALPHA_4 = 0.0,
                      BETA_4  = 1.0;
*/

static constexpr ldbl A = 0.0,
                      B = 0.0,
                      C = 0.0;
static constexpr ldbl L_1 = 1.0,
                      L_2 = 1.0;
static constexpr ldbl ALPHA_1 = 0.0,
                      BETA_1  = 1.0,
                      ALPHA_2 = 0.0,
                      BETA_2  = 1.0,
                      ALPHA_3 = 0.0,
                      BETA_3  = 1.0,
                      ALPHA_4 = 0.0,
                      BETA_4  = 1.0;

int main(int argc, const char *argv[])
{
    try
    {
        ldbl epsilon, y;
        std::size_t n_1, n_2;
        std::cin >> n_1 >> n_2 >> epsilon >> y;

        const ldbl h_1 = L_1 / n_1, h_2 = L_2 / n_2;
        std::vector<ldbl> liebmann_fdm_error,
            seidel_fdm_error;
        const ublas::matrix<ldbl>
            liebmann_fdm_u = liebmann_fdm<ldbl>(A, B, C, f_x_t,
                L_1, n_1, L_2, n_2, ALPHA_1, BETA_1, ALPHA_2, BETA_2,
                ALPHA_3, BETA_3, ALPHA_4, BETA_4,
                phi_1_y, phi_2_y, phi_3_x, phi_4_x,
                [&] (const ublas::matrix<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size1 = u_k.size1(), size2 = u_k.size2();
                    for (std::size_t i = 0; i < size1; ++i)
                    {
                        for (std::size_t j = 0; j < size2; ++j)
                        {
                            error = std::max(error, std::abs(u_k(i, j) -
                                u_exact(i * h_1, j * h_2)));
                        }
                    }
                    seidel_fdm_error.push_back(error);
                }, epsilon),
            seidel_fdm_u = seidel_fdm<ldbl>(A, B, C, f_x_t,
                L_1, n_1, L_2, n_2, ALPHA_1, BETA_1, ALPHA_2, BETA_2,
                ALPHA_3, BETA_3, ALPHA_4, BETA_4,
                phi_1_y, phi_2_y, phi_3_x, phi_4_x,
                [&] (const ublas::matrix<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size1 = u_k.size1(), size2 = u_k.size2();
                    for (std::size_t i = 0; i < size1; ++i)
                    {
                        for (std::size_t j = 0; j < size2; ++j)
                        {
                            error = std::max(error, std::abs(u_k(i, j) -
                                u_exact(i * h_1, j * h_2)));
                        }
                    }
                    liebmann_fdm_error.push_back(error);
                }, epsilon);
/*
        std::cout << std::round(y / L_2 * n_2) << '\n';
        ublas::matrix<ldbl> u_k_exact(n_1 + 1, n_2 + 1);
        for (std::size_t i = 0; i <= n_1; ++i)
            for (std::size_t j = 0; j <= n_2; ++j)
                u_k_exact(i, j) = u_exact(i * h_1, j * h_2);
        std::cout << "u_k_exact = ";
        print_m(u_k_exact);
*/
        const std::size_t j = std::round(y / L_2 * n_2);
        for (std::size_t i = 0; i < liebmann_fdm_u.size1(); ++i)
        {
            std::cout << std::fixed << std::setprecision(18) << i * h_1 << ',' <<
                u_exact(i * h_1, j * h_2) << ',' << liebmann_fdm_u(i, j) << ',' <<
                seidel_fdm_u(i, j) << '\n';
        }
        if (argc > 1)
        {
            std::fstream error_stream(argv[1], std::ios_base::out | std::ios_base::trunc);
            for (std::size_t i = 0; i < liebmann_fdm_error.size(); ++i)
            {
                error_stream << std::fixed << std::setprecision(18) << i << ',' <<
                    liebmann_fdm_error[i] << ',' << seidel_fdm_error[i] << '\n';
            }
        }
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
/*
static constexpr ldbl f_x_t(const ldbl &, const ldbl &) noexcept
{
    return 0.0;
}

static constexpr ldbl phi_1_y(const ldbl &y) noexcept
{
    return std::exp(-y) * std::cos(y);
}

static constexpr ldbl phi_2_y(const ldbl &) noexcept
{
    return 0.0;
}

static constexpr ldbl phi_3_x(const ldbl &x) noexcept
{
    return std::exp(-x) * std::cos(x);
}

static constexpr ldbl phi_4_x(const ldbl &) noexcept
{
    return 0.0;
}

static constexpr ldbl u_exact(const ldbl &x, const ldbl &y) noexcept
{
    return std::exp(-x - y) * std::cos(x) * std::cos(y);
}
*/
static constexpr ldbl f_x_t(const ldbl &, const ldbl &) noexcept
{
    return 0.0;
}

static constexpr ldbl phi_1_y(const ldbl &y) noexcept
{
    return y;
}

static constexpr ldbl phi_2_y(const ldbl &y) noexcept
{
    return 1.0 + y;
}

static constexpr ldbl phi_3_x(const ldbl &x) noexcept
{
    return x;
}

static constexpr ldbl phi_4_x(const ldbl &x) noexcept
{
    return 1.0 + x;
}

static constexpr ldbl u_exact(const ldbl &x, const ldbl &y) noexcept
{
    return x + y;
}
