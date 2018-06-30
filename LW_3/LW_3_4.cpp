#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "derivation.hpp"

typedef double dbl;

template <typename T>
void derivation(std::ostream &, std::vector<T> &x, std::vector<T> &f_of_x, const T);

template <typename T>
void in_x_and_f_of_x(std::istream &, std::vector<T> &, std::vector<T> &);

int main()
{
    dbl x_star;
    std::cin >> x_star;
    std::vector<dbl> x, f_of_x;
    in_x_and_f_of_x(std::cin, x, f_of_x);
    try
    {
        derivation(std::cout, x, f_of_x, x_star);
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }
    std::cout << std::endl;

    return 0;
}

template <typename T>
void in_x_and_f_of_x(std::istream &is, std::vector<T> &x, std::vector<T> &f_of_x)
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

template <typename T>
void derivation(std::ostream &os, std::vector<T> &x, std::vector<T> &f_of_x, const T x_star)
{
    const std::size_t size = x.size();
    std::sort(x.begin(), x.end());
    const typename std::vector<T>::const_iterator begin = x.cbegin(), end = x.cend(),
        it = std::lower_bound(begin, end, x_star);
    if (size < 2 || x_star < *begin || *(end - 1) < x_star)
    {
        throw std::domain_error("Out x out of range");
    }
    os << "x = " << x_star << '\n';
    T der;
    bool enough_in_left = (it != begin), enough_in_right = (it + 1 != end);
    if (x_star == *it)
    {
        if (enough_in_left)
        {
            der = df_by_dx_order_1(x, f_of_x, it - begin - 1);
            os << "left  side derivation: f'(x - 0) = " << der << '\n';
        }
        if (enough_in_right)
        {
            der = df_by_dx_order_1(x, f_of_x, it - begin);
            os << "right side derivation: f'(x + 0) = " << der << '\n';
        }
    }
    else
    {
        der = df_by_dx_order_1(x, f_of_x, it - begin - 1);
        os << "derivation: f'(x - 0) = " << der << '\n';
    }
    if (enough_in_left && enough_in_right)
    {
        der = df_by_dx_order_2(x, f_of_x, x_star, it - begin - 1);
        os << "derivation (2-precision):  f'(x) = " << der << '\n';
        der = d2f_by_dx2(x, f_of_x, it - begin - 1);
        os << "second derivation:         f'(x) = " << der << '\n';
    }
}
