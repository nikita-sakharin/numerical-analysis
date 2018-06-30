#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "least_squares_method.hpp"

enum Order : std::size_t
{
    ONE = 1,
    TWO,
    THREE
};

template <typename T, std::size_t ORDER>
using Polynom = std::array<T, ORDER + 1>;

typedef double dbl;

template <typename T>
void in_x_and_f_of_x(std::istream &, ublas::vector<T> &, ublas::vector<T> &);

template<typename T, std::size_t ORDER,
    typename = std::enable_if_t<ORDER>>
std::ostream &polynom_out(std::ostream &, const Polynom<T, ORDER> &);

int main()
{
    ublas::vector<dbl> x, f_of_x;
    in_x_and_f_of_x(std::cin, x, f_of_x);

    try
    {
        Polynom<dbl, ONE> p1 = std::move(least_squares_method<dbl, ONE>(x, f_of_x));
        polynom_out<dbl, ONE>(std::cout, p1);
        Polynom<dbl, TWO> p2 = std::move(least_squares_method<dbl, TWO>(x, f_of_x));
        polynom_out<dbl, TWO>(std::cout, p2);

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

template<typename T, std::size_t ORDER,
    typename = std::enable_if_t<ORDER>>
std::ostream &polynom_out(std::ostream &os, const Polynom<T, ORDER> &polynom)
{
    std::cout << polynom[0];
    for (std::size_t i = 0; i < ORDER; ++i)
    {
        std::cout << " + " << polynom[i] << "* x^" << i + 1;
    }
    std::cout << '\n';

    return os;
}
