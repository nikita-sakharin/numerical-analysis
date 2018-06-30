#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "jacobi_eigenvalue_algorithm.hpp"

int main()
{
    try
    {
        Matrix<long double> a;
        std::cin >> a;

        std::pair<Matrix<long double>, Matrix<long double>>
            result = jacobi_eigenvalue_algorithm(a, 0.01);
        const typename Matrix<long double>::size_type m = a.size1();
        std::cout << "eigenvalue:\n";
        for (typename Matrix<long double>::size_type i = 0u; i < m; ++i)
        {
            std::cout << result.first(i, i) << '\n';
        }
        std::cout << '\n';
        std::cout << result.second << '\n';
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
