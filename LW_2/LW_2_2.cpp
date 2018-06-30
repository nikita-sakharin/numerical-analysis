#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <exception>

#include "root_finding_system.hpp"

using namespace boost::numeric;

typedef double ldbl;

static constexpr const ldbl A = 3.0;

ublas::vector<ldbl> f(const ublas::vector<ldbl> &);
ublas::matrix<ldbl> jacobi(const ublas::vector<ldbl> &);
ublas::vector<ldbl> g(const ublas::vector<ldbl> &);

int main()
{
    try
    {
        ublas::vector<ldbl> x0(2);
        x0(0) = x0(1) = 1.0;
        ublas::vector<ldbl> root;

        root = newton_method_system<ldbl>(f, jacobi, x0, 0.00001);
        std::cout << "x* = " << root << '\n';
        std::cout << "f(x*) = " << f(root) << '\n';

        root = successive_iteration_system<ldbl>(g, x0, 0.00001);
        std::cout << "x* = " << root << '\n';
        std::cout << "f(x*) = " << f(root) << '\n';
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

ublas::vector<ldbl> f(const ublas::vector<ldbl> &x)
{
    ublas::vector<ldbl> f_of_x(2);
    f_of_x(0) = x(0) - std::cos(x(1)) - 1.0;
    f_of_x(1) = x(1) - std::log10(x(0) + 1.0) - A;

    return f_of_x;
}

ublas::matrix<ldbl> jacobi(const ublas::vector<ldbl> &x)
{
    ublas::matrix<ldbl> jacobi_of_x(2, 2);
    jacobi_of_x(0, 0) = 1.0;
    jacobi_of_x(0, 1) = std::sin(x(1));
    jacobi_of_x(1, 0) = 1.0 / (x(0) + 1.0) / log(10.0);
    jacobi_of_x(1, 1) = 1.0;

    return jacobi_of_x;
}

ublas::vector<ldbl> g(const ublas::vector<ldbl> &x)
{
    ublas::vector<ldbl> g_of_x(2);
    g_of_x(0) = std::cos(x(1)) + 1.0;
    g_of_x(1) = std::log10(x(0) + 1.0) + A;

    return g_of_x;
}
