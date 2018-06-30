#ifndef __INTERPOLATION_METHOD_HPP__
#define __INTERPOLATION_METHOD_HPP__

#include <cstddef>

#include <type_traits>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric;

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T lagrange_method(const ublas::vector<T> &, const ublas::vector<T> &, const T,
    ublas::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T newton_method(const ublas::vector<T> &, const ublas::vector<T> &, const T,
    ublas::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T divided_diff(const std::size_t, const std::size_t,
    const ublas::vector<T> &, const ublas::vector<T> &);

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T lagrange_method(const ublas::vector<T> &x_i, const ublas::vector<T> &f_of_x_i,
    const T x, ublas::vector<T> &polynom)
{
    using namespace std;
    const size_t size = x_i.size();
    polynom.resize(size);
    T f_of_x = 0., w = 1.;
    for (size_t i = 0; i < size; ++i)
    {
        w *= (x - x_i[i]);
        polynom[i] = f_of_x_i[i];
    }
    for (size_t i = 0; i < size; ++i)
    {
        T prod = w;
        for (size_t j = 0; j < size; ++j)
        {
            if (j == i)
            {
                prod /= (x - x_i[i]);
            }
            else
            {
                const T divider = x_i[i] - x_i[j];
                polynom[i] /= divider;
                prod /= divider;
            }
        }
        f_of_x += f_of_x_i[i] * prod;
    }

    return f_of_x;
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T newton_method(const ublas::vector<T> &x_i, const ublas::vector<T> &f_of_x_i,
    const T x, ublas::vector<T> &polynom)
{
    using namespace std;
    const size_t size = x_i.size();
    polynom.resize(size);
    T prod = x - x_i[0], f_of_x = polynom[0] = f_of_x_i[0];
    for (size_t i = 1; i < size; ++i)
    {
        const T diff = divided_diff(0, i, x_i, f_of_x_i);
        polynom[i] = diff;
        f_of_x += prod * diff;
        prod *= x - x_i[i];
    }

    return f_of_x;
}

template<typename T,
    typename = std::enable_if_t<std::is_floating_point<T>::value>>
T divided_diff(const std::size_t i, const std::size_t j,
    const ublas::vector<T> &x_i, const ublas::vector<T> &f_of_x_i)
{
    if (i > j)
    {
        return divided_diff(j, i, x_i, f_of_x_i);
    }
    else if (i == j)
    {
        return f_of_x_i[i];
    }

    return (divided_diff(i, j - 1, x_i, f_of_x_i) -
        divided_diff(i + 1, j, x_i, f_of_x_i)) / (x_i[i] - x_i[j]);
}

#endif
