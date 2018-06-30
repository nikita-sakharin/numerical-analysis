#ifndef __ROOT_FINDING_SYSTEM_HPP__
#define __ROOT_FINDING_SYSTEM_HPP__

#include <cmath>

#include "../generic/lup.hpp"

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
ublas::vector<T> newton_method_system(ublas::vector<T> (const ublas::vector<T> &),
    ublas::matrix<T> (const ublas::vector<T> &),
    const ublas::vector<T> &, const T);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
ublas::vector<T> successive_iteration_system(ublas::vector<T> f(const ublas::vector<T> &),
    const ublas::vector<T> &x0, const T epsilon);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
ublas::vector<T> newton_method_system(ublas::vector<T> f(const ublas::vector<T> &),
    ublas::matrix<T> jacobi(const ublas::vector<T> &),
    const ublas::vector<T> &x0, const T epsilon)
{
    ublas::vector<T> x_k, x_k_plus_1 = x0, delta_x;
    do
    {
        x_k = x_k_plus_1;
        delta_x = LUP<T>(jacobi(x_k)).solution(-f(x_k));
        x_k_plus_1 = x_k + delta_x;
    } while (norm_2(x_k_plus_1 - x_k) >= epsilon);

    return x_k_plus_1;
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
ublas::vector<T> successive_iteration_system(ublas::vector<T> f(const ublas::vector<T> &),
    const ublas::vector<T> &x0, const T epsilon)
{
    ublas::vector<T> x_k, x_k_plus_1 = x0;
    do
    {
        x_k = x_k_plus_1;
        x_k_plus_1 = f(x_k);
    } while (norm_2(x_k_plus_1 - x_k) >= epsilon);

    return x_k_plus_1;
}

#endif
