#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "LUP.hpp"

int main()
{
    try
    {
        Matrix<double> m;
        std::cin >> m;
        LUP<double> lup(m);
        std::cout << lup.invertible() << '\n';

    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
