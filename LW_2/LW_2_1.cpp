#include <cmath>

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "root_finding.hpp"

typedef double dbl;

dbl f(dbl);
dbl f_der(dbl);
dbl g(dbl);

int main()
{
    std::cout << newton_method<dbl>(f, f_der, 3.0, 0.001) << '\n';
    std::cout << successive_iteration<dbl>(g, 3.0, 0.001) << '\n';

    return 0;
}

dbl f(dbl x)
{
    return std::exp(x) - 2.0 * x - 2.0;
}

dbl f_der(dbl x)
{
    return std::exp(x) - 2.0;
}

dbl g(dbl x)
{
    return std::log(2.0 * x + 2.0);
}
