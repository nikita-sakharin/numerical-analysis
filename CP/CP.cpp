#include <cmath>

#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "lanczos_algorithm.hpp"

typedef long double ldbl;

int main()
{

    try
    {
        ublas::matrix<ldbl> a_in;
        ublas::vector<ldbl> b_in;
        std::cin >> a_in >> b_in;
        SparseMatrix<ldbl> a = a_in;
        SparseVector<ldbl> b = b_in;
        SparseVector<ldbl> eigenvalue = lanczos_algorithm(a, b);
        const std::size_t size = eigenvalue.size();
        std::cout << "eigenvalue : \n";
        for (std::size_t i = 0; i < size; ++i)
        {
            std::cout << eigenvalue(i) << '\n';
        }
    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
