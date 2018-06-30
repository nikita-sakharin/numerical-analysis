#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "lup.hpp"

int main()
{
    try
    {
        Matrix<long double> m, b;
        std::cin >> m >> b;
        LUP<long double> lup(m);
        std::cout << "L:\n" << lup.get_l() << '\n';
        std::cout << "U:\n" << lup.get_u() << '\n';
        std::cout << "P:\n" << lup.get_p() << '\n';
        std::cout << "Solution:\n" << lup.solution(b) << '\n';
        std::cout << "Invertible:\n" << lup.invertible() << '\n';
        std::cout << "A * A ^ (-1):\n" << lup.invertible() * m << '\n';

    } catch (const std::exception &ex)
    {
        std::cout << ex.what() << '\n';
    }

    return 0;
}
