#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "qr.hpp"

int main()
{
    try
    {
        Matrix<long double> a;
        std::cin >> a;
        Matrix<std::complex<long double>> result = qr(a, 0.01);
        const typename Matrix<long double>::size_type m = a.size1();
        std::cout << "eigenvalue:\n";
        for (typename Matrix<long double>::size_type i = 0u; i < m; ++i)
        {
            std::cout << result(i) << '\n';
        }
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
