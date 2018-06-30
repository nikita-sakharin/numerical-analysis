#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "iterative_method.hpp"

int main()
{
    try
    {
        Matrix<long double> a, b, x;
        std::cin >> a >> b;
        iteration_count k, k_real;

        x = successive_iteration<long double, order_infinity>(a, b, 0.01, k, k_real);
        std::cout << "k_theoretical = " << k << '\n';
        std::cout << "k_real = " << k_real << '\n';
        std::cout << x << '\n';

        x = gauss_seidel<long double, order_infinity>(a, b, 0.01, k, k_real);
        std::cout << "k_theoretical = " << k << '\n';
        std::cout << "k_real = " << k_real << '\n';
        std::cout << x << '\n';
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
