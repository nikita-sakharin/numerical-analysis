#include <cmath>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "interpolation_method.hpp"

typedef double dbl;

dbl f(dbl);

template <typename T>
void test(std::ostream &, std::istream &, const T);

template <typename T>
void in_x(std::istream &, ublas::vector<T> &);

template <typename T>
void compute_f_of_x(const ublas::vector<T> &, ublas::vector<T> &, T f(T));

template <typename T>
std::ostream &lagrange_polynom_out(std::ostream &, const ublas::vector<T> &);

template <typename T>
std::ostream &newton_polynom_out(std::ostream &, const ublas::vector<T> &,
    const ublas::vector<T> &);

int main()
{
    dbl x_star;
    std::cin >> x_star;
    const dbl f_of_x_star = f(x_star);
    std::cout << "f(x*): " << f_of_x_star << '\n' << std::endl;

    test(std::cout, std::cin, x_star);
    test(std::cout, std::cin, x_star);

    return 0;
}

template <typename T>
void test(std::ostream &os, std::istream &is, const T x_star)
{
    T f_value;
    ublas::vector<dbl> x, f_of_x, polynom;
    in_x(is, x);
    compute_f_of_x(x, f_of_x, f);

    f_value = lagrange_method<dbl>(x, f_of_x, x_star, polynom);
    lagrange_polynom_out(os, polynom);
    os << "lagrange_method: P(" << x_star << ") = " << f_value << '\n' << std::endl;

    f_value = newton_method<dbl>(x, f_of_x, x_star, polynom);
    newton_polynom_out(os, x, polynom);
    os << "newton_method:   P(" << x_star << ") = " << f_value << '\n' << std::endl;
}

dbl f(dbl x)
{
    return exp(x);
}

template <typename T>
void in_x(std::istream &is, ublas::vector<T> &x)
{
    std::size_t n;
    is >> n;
    x.resize(n);
    for (T &x_i : x)
    {
        is >> x_i;
    }
}

template <typename T>
void compute_f_of_x(const ublas::vector<T> &x, ublas::vector<T> &f_of_x, T f(T))
{
    f_of_x.resize(x.size());
    std::transform(x.cbegin(), x.cend(), f_of_x.begin(), f);
}

template <typename T>
std::ostream &lagrange_polynom_out(std::ostream &os, const ublas::vector<T> &p)
{
    const std::size_t size = p.size();
    os << "P(x) =";
    for (std::size_t i = 0; i < size; ++i)
    {
        os << ' ' << p[i];
        os << std::showpos;
        for (std::size_t j = 0; j < size; ++j)
        {
            if (i != j)
            {
                os << "(x " << p[j] << ")";
            }
        }
    }
    os << std::noshowpos << '\n';

    return os;
}

template <typename T>
std::ostream &newton_polynom_out(std::ostream &os, const ublas::vector<T> &x,
    const ublas::vector<T> &p)
{
    static const std::string sign = "+-";
    const std::size_t size = p.size();
    std::string prod;
    os << "P(x) =";
    for (std::size_t i = 0; i < size; ++i)
    {
        os << ' ' << p[i] << prod;
        os << std::showpos;
        prod += "(x " + (sign[x[i] >= 0.0] +  std::to_string(fabs(x[i]))) + ")";
    }
    os << std::noshowpos << '\n';

    return os;
}
