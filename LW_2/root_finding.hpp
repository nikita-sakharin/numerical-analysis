#ifndef __ROOT_FINDING_HPP__
#define __ROOT_FINDING_HPP__

#include <cmath>

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T newton_method(T (T), T (T), T, T);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T successive_iteration(T (T), T, T);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T newton_method(T f(T), T f_der(T), const T x0, const T epsilon)
{
    T x_k, x_k_plus_1 = x0;
    do
    {
        x_k = x_k_plus_1;
        x_k_plus_1 = x_k - f(x_k) / f_der(x_k);
    } while (std::fabs(x_k_plus_1 - x_k) >= epsilon);

    return x_k_plus_1;
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T successive_iteration(T f(T), const T x0, const T epsilon)
{
    T x_k, x_k_plus_1 = x0;
    do
    {
        x_k = x_k_plus_1;
        x_k_plus_1 = f(x_k);
    } while (std::fabs(x_k_plus_1 - x_k) >= epsilon);

    return x_k_plus_1;
}

#endif
