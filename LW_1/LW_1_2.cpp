#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "thomas_algorithm.hpp"

template <typename T>
void get_a_b_c_d(std::istream &, Matrix<T> &, Matrix<T> &, Matrix<T> &,
    Matrix<T> &);

int main()
{
    try
    {
        Matrix<long double> a, b, c, d;
        get_a_b_c_d(std::cin, a, b, c, d);
        std::cout << "Solution:\n" << thomas_algorithm(a, b, c, d) << '\n';
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}

template <typename T>
void get_a_b_c_d(std::istream &is, Matrix<T> &a, Matrix<T> &b, Matrix<T> &c,
    Matrix<T> &d)
{
    typename Matrix<T>::size_type n;
    is >> n;
    a.resize(1u, n);
    b.resize(1u, n);
    c.resize(1u, n);
    d.resize(1u, n);

    for (typename Matrix<T>::size_type i = 0u; i < n; ++i)
    {
        if (i)
        {
            is >> a(i);
        }
        else
        {
            a(i) = 0.0;
        }
        is >> b(i);
        if (i != n - 1u)
        {
            is >> c(i);
        }
        else
        {
            c(i) = 0.0;
        }

        is >> d(i);
    }
}
