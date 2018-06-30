#ifndef __LEAST_SQUARES_METHOD_HPP__
#define __LEAST_SQUARES_METHOD_HPP__

#include <cstddef>

#include <array>
#include <type_traits>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../generic/lup.hpp"

using namespace boost::numeric;

template <typename T, std::size_t ORDER>
using Polynom = std::array<T, ORDER + 1>;

template<typename T, std::size_t ORDER,
    typename = std::enable_if_t<std::is_floating_point<T>::value && ORDER>>
Polynom<T, ORDER> least_squares_method(const ublas::vector<T> &,
    const ublas::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T element_prod_sum(ublas::vector<T> &, const ublas::vector<T> &);

template<typename T, std::size_t ORDER,
    typename = std::enable_if_t<std::is_floating_point<T>::value && ORDER>>
Polynom<T, ORDER> least_squares_method(const ublas::vector<T> &x,
    const ublas::vector<T> &f_of_x)
{
    static constexpr const std::size_t ORDER_PLUS_1 = ORDER + 1u;
    const std::size_t size = x.size();

    if (!size || size != f_of_x.size())
    {
        throw std::domain_error("!x.size() || x.size() != f_of_x.size()");
    }

    ublas::matrix<T> a(ORDER_PLUS_1, ORDER_PLUS_1);
    ublas::vector<T> b(ORDER_PLUS_1), curr_x(size, 1.0), curr_y_x = f_of_x;
    T sum_x = static_cast<T>(size), sum_y_x = ublas::sum(f_of_x);

    const std::ptrdiff_t diag_cnt =
        static_cast<std::ptrdiff_t>(ORDER_PLUS_1 * 2u - 1u);
    for (std::ptrdiff_t k = 0; k < diag_cnt; ++k)
    {
        std::ptrdiff_t begin = k, end = -1;
        if (k < static_cast<std::ptrdiff_t>(ORDER_PLUS_1))
        {
            b(k) = sum_y_x;
        }
        else
        {
            begin = ORDER_PLUS_1 - 1;
            end = k - ORDER_PLUS_1;
        }
        for (std::ptrdiff_t j = begin; j > end; --j)
        {
            a(k - j, j) = sum_x;
        }
        sum_x = element_prod_sum(curr_x, x);
        sum_y_x = element_prod_sum(curr_y_x, x);
    }

    LUP<T> lup(a);
    ublas::vector<T> x_star = LUP<T>(a).solution(b);

    Polynom<T, ORDER> result;
    std::copy_n(x_star.cbegin(), ORDER_PLUS_1, result.begin());

    return result;
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T element_prod_sum(ublas::vector<T> &lhs, const ublas::vector<T> &rhs)
{
    const std::size_t size = lhs.size();
    T sum = 0.0;
    for (std::size_t i = 0; i < size; ++i)
    {
        lhs[i] *= rhs[i];
        sum += lhs[i];
    }
    return sum;
}

#endif
