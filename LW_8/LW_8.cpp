#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../generic/header.hpp"
#include "splitting_methods.hpp"

static constexpr ldbl f_x_t(const ldbl &, const ldbl &) noexcept;
static constexpr ldbl phi_1_y(const ldbl &) noexcept;
static constexpr ldbl phi_2_y(const ldbl &) noexcept;
static constexpr ldbl phi_3_x(const ldbl &) noexcept;
static constexpr ldbl phi_4_x(const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &) noexcept;

static constexpr ldbl QUIET_NAN_LDBL = std::numeric_limits<ldbl>::quiet_NaN();

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

int main(int argc, const char *argv[])
{
    try
    {
        ldbl epsilon, y, omega;
        std::size_t n_1, n_2;
        std::cin >> n_1 >> n_2 >> epsilon >> y >> omega;

        const ldbl h_1 = L_1 / n_1, h_2 = L_2 / n_2;
        std::vector<ldbl> successive_error, seidel_error;
        const ublas::matrix<ldbl>
            successive_fdm_u = successive_fdm<ldbl>(A, B, C, f_x_t,
                L_1, n_1, L_2, n_2, ALPHA_1, BETA_1, ALPHA_2, BETA_2,
                ALPHA_3, BETA_3, ALPHA_4, BETA_4,
                phi_1_y, phi_2_y, phi_3_x, phi_4_x,
                [&] (const ublas::matrix<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size1 = u_k.size1() - 1, size2 = u_k.size2() - 1;
                    for (std::size_t i = 1; i < size1; ++i)
                    {
                        for (std::size_t j = 1; j < size2; ++j)
                        {
                            error = std::max(error, std::abs(u_k(i, j) -
                                u_exact(i * h_1, j * h_2)));
                        }
                    }
                    successive_error.push_back(error);
                }, epsilon),
            seidel_fdm_u = seidel_fdm<ldbl>(A, B, C, f_x_t,
                L_1, n_1, L_2, n_2, ALPHA_1, BETA_1, ALPHA_2, BETA_2,
                ALPHA_3, BETA_3, ALPHA_4, BETA_4,
                phi_1_y, phi_2_y, phi_3_x, phi_4_x, omega,
                [&] (const ublas::matrix<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const std::size_t size1 = u_k.size1() - 1, size2 = u_k.size2() - 1;
                    for (std::size_t i = 1; i < size1; ++i)
                    {
                        for (std::size_t j = 1; j < size2; ++j)
                        {
                            error = std::max(error, std::abs(u_k(i, j) -
                                u_exact(i * h_1, j * h_2)));
                        }
                    }
                    seidel_error.push_back(error);
                }, epsilon);
        const std::size_t j = std::round(y / L_2 * n_2);
        for (std::size_t i = 0; i < successive_fdm_u.size1(); ++i)
        {
            std::cout << std::fixed << std::setprecision(18) << i * h_1 << ',' <<
                u_exact(i * h_1, j * h_2) << ',' << successive_fdm_u(i, j) << ',' <<
                seidel_fdm_u(i, j) << '\n';
        }
        if (argc > 1)
        {
            std::fstream error_stream(argv[1], std::ios_base::out | std::ios_base::trunc);
            for (std::size_t i = 0;
                i < std::max(successive_error.size(), seidel_error.size()); ++i)
            {
                error_stream << std::fixed << std::setprecision(18) << i;
                error_stream << ','
                    << (i < successive_error.size() ? successive_error[i] : QUIET_NAN_LDBL);
                error_stream << ','
                    << (i < seidel_error.size() ? seidel_error[i] : QUIET_NAN_LDBL);
                error_stream << '\n';
            }
        }
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

static constexpr ldbl f_x_t(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &x, const ldbl &y, const ldbl &t) noexcept
{
    return 0.0;
}

static constexpr ldbl phi_1_x_t(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &x, const ldbl &t) noexcept
{
    return ;
}

static constexpr ldbl phi_2_x_t(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &x, const ldbl &t) noexcept
{
    return ;
}

static constexpr ldbl phi_3_y_t(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &y, const ldbl &t) noexcept
{
    return ;
}

static constexpr ldbl phi_4_y_t(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &y, const ldbl &t) noexcept
{
    return -std::sin(x) * std::sin(mu * t)
}

static constexpr ldbl psi_x_y(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &x, const ldbl &y) noexcept
{
    return 0.0;
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &b, const ldbl &mu,
    const ldbl &x, const ldbl &y, const ldbl &t) noexcept
{
    return std::sin(x) * std::sin(y) * std::sin(mu * t);
}
