#ifndef __DERIVATION_HPP__
#define __DERIVATION_HPP__

#include <cstddef>

#include <type_traits>
#include <vector>

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T df_by_dx_order_1(const std::vector<T> &x, const std::vector<T> &,
    std::size_t);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T df_by_dx_order_2(const std::vector<T> &x, const std::vector<T> &, const T,
    std::size_t);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T d2f_by_dx2(const std::vector<T> &x, const std::vector<T> &, T,
    std::size_t);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T df_by_dx_order_1(const std::vector<T> &x, const std::vector<T> &f_of_x,
    const std::size_t i)
{
    const std::size_t size = x.size(), i_plus_1 = i + 1;
    if (size < 2 || 1 > i_plus_1 || i_plus_1 >= size || size != f_of_x.size())
    {
        throw std::domain_error("Can't be calculated");
    }

    return (f_of_x[i_plus_1] - f_of_x[i]) / (x[i_plus_1] - x[i]);
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T df_by_dx_order_2(const std::vector<T> &x, const std::vector<T> &f_of_x,
    const T x_star, const std::size_t i)
{
    const std::size_t size = x.size(), i_plus_1 = i + 1, i_plus_2 = i + 2;
    if (size < 3 || 2 > i_plus_2 || i_plus_2 >= size || size != f_of_x.size())
    {
        throw std::domain_error("Can't be calculated");
    }

    const T term = (f_of_x[i_plus_1] - f_of_x[i]) / (x[i_plus_1] - x[i]),
        first_term = (f_of_x[i_plus_2] - f_of_x[i_plus_1]) / (x[i_plus_2] - x[i_plus_1]),
        second_term = (f_of_x[i_plus_1] - f_of_x[i]) / (x[i_plus_1] - x[i]);
    return term + (first_term - second_term) / (x[i_plus_2] - x[i]) *
        (2.0 * x_star - x[i] - x[i_plus_1]);
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T d2f_by_dx2(const std::vector<T> &x, const std::vector<T> &f_of_x,
    const std::size_t i)
{
    const std::size_t size = x.size(), i_plus_1 = i + 1, i_plus_2 = i + 2;
    if (size < 3 || 2 > i_plus_2 || i_plus_2 >= size || size != f_of_x.size())
    {
        throw std::domain_error("Can't be calculated");
    }

    const T first_term = (f_of_x[i_plus_2] - f_of_x[i_plus_1]) / (x[i_plus_2] - x[i_plus_1]),
        second_term = (f_of_x[i_plus_1] - f_of_x[i]) / (x[i_plus_1] - x[i]);
    return 2.0 * (first_term - second_term) / (x[i_plus_2] - x[i]);
}

#endif
