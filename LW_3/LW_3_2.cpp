#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "cubic_spline.hpp"

typedef double dbl;

template <typename T>
void in_x_and_f_of_x(std::istream &, ublas::vector<T> &, ublas::vector<T> &);

int main()
{
    dbl x_star;
    std::cin >> x_star;
    ublas::vector<dbl> x, f_of_x;

    in_x_and_f_of_x(std::cin, x, f_of_x);
    try
    {
        std::cout << "cubic_spline in x = " << x_star
            << ", f(x) = " << cubic_spline(std::cout, x, f_of_x, x_star);
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }
    std::cout << std::endl;

    return 0;
}

template <typename T>
void in_x_and_f_of_x(std::istream &is, ublas::vector<T> &x, ublas::vector<T> &f_of_x)
{
    std::size_t n;
    is >> n;
    x.resize(n);
    f_of_x.resize(n);

    for (T &x_i : x)
    {
        is >> x_i;
    }
    for (T &f_of_x_i : f_of_x)
    {
        is >> f_of_x_i;
    }
}
